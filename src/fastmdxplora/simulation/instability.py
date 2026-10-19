"""Observe finite forces in the running Context; retain evidence only on failure.

The observer neither edits the Context nor retries integration. One finite
observation is held in memory. A force rank is an observation, not a cause.
"""
from __future__ import annotations

import hashlib
import io
import json
import logging
import re
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)
MAX_ATOMS = 200_000
MAX_SNAPSHOT = 32 * 1024 * 1024
RANKS = 12
NOTICE = ("Last recorded finite observation, not necessarily the step before failure. "
          "Large forces do not establish the cause. Forces include the active "
          "integration force groups and any restraints or bias in those groups.")


def numeric_failure(error: BaseException) -> bool:
    """Keep numerical failures distinct from an I/O error or a requested stop."""
    seen = set()
    while error is not None and id(error) not in seen:
        seen.add(id(error))
        if re.search(r"\b(?:nan|inf|infinite|non-finite)\b", str(error), re.I):
            return True
        error = error.__cause__ or error.__context__
    return False


class FiniteForceObserver:
    """OpenMM's reporter protocol, compatible with its tuple-based releases.

    Captures full-system positions and the largest finite force norms from
    the *same* State. Disk and snapshot formatting work happen only on a
    numerical failure. Capture failures are best effort and never escape.
    """

    def __init__(self, simulation: Any, output: Path, interval: int, *, supported: bool = True):
        self.output = Path(output)
        self.interval = max(1, int(interval))
        self.stage = "initial"
        self.last = None
        self.problem = "No finite force observation was retained."
        self.model = ""
        self.groups = None
        self.enabled = False
        self.config_digest = ""
        self.publish = True
        try:
            old = self.output / "instability.json"
            if old.exists() or old.is_symlink():
                # A deliberately rerun phase must not show a previous run's
                # observation as its own, or destroy the earlier evidence.
                if old.is_symlink() or old.stat().st_size > 100_000:
                    from uuid import uuid4

                    old.replace(self.output / f"instability-unreadable-{uuid4().hex}.json")
                else:
                    raw = old.read_bytes()
                    digest = hashlib.sha256(raw).hexdigest()
                    old.replace(self.output / f"instability-record-{digest}.json")
        except OSError:
            self.publish = False
        try:
            config = self.output.parent / "resolved_config.yml"
            if config.stat().st_size <= 1_000_000:
                self.config_digest = hashlib.sha256(config.read_bytes()).hexdigest()
        except OSError:
            pass
        try:
            self.atoms = list(simulation.topology.atoms())
            self.enabled = supported and 0 < len(self.atoms) <= MAX_ATOMS
            if not supported:
                self.problem = "Finite force capture is unavailable for configured enhanced sampling: extra force observations have not been verified to preserve its dynamics."
            elif not self.enabled:
                self.problem = "The full system exceeds the force-evidence atom limit."
        except Exception:
            self.problem = "The full-system topology could not be read."

    def begin(self, simulation: Any, stage: str) -> None:
        self.stage = stage
        self.last = None  # a previous stage's force model is not this one
        if not self.enabled:
            return
        try:
            from openmm import XmlSerializer

            self.model = hashlib.sha256(XmlSerializer.serialize(simulation.system).encode()).hexdigest()
            self.groups = int(simulation.integrator.getIntegrationForceGroups())
            state = simulation.context.getState(getPositions=True, getForces=True,
                getEnergy=True, getParameters=True, enforcePeriodicBox=False, groups=self.groups)
            self.report(simulation, state)
        except Exception:
            self.problem = "A finite force observation could not be obtained at this stage's start."

    def describeNextReport(self, simulation: Any):
        return (self.interval - simulation.currentStep % self.interval,
                True, False, True, True, False)

    def report(self, simulation: Any, state: Any) -> None:
        if not self.enabled:
            return
        try:
            import numpy as np
            from openmm import unit

            xyz = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
            forces = np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer))
            energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
            time = float(state.getTime().value_in_unit(unit.picosecond))
            box = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
            if xyz.shape != (len(self.atoms), 3) or forces.shape != xyz.shape:
                self.problem = "Force and full-system position counts did not match."
                return
            if not (np.isfinite(xyz).all() and np.isfinite(forces).all()
                    and np.isfinite(box).all() and np.isfinite([energy, time]).all()):
                return  # a bad observation must never replace the last finite one
            norms = np.hypot(np.hypot(forces[:, 0], forces[:, 1]), forces[:, 2])
            if not np.isfinite(norms).all():
                return
            count = min(RANKS, len(norms))
            threshold = np.partition(norms, len(norms) - count)[len(norms) - count]
            above = np.flatnonzero(norms > threshold)
            tied = np.flatnonzero(norms == threshold)[:count - len(above)]
            indices = np.concatenate((above, tied))
            indices = indices[np.lexsort((indices, -norms[indices]))]
            parameters = {str(key): float(value) for key, value in simulation.context.getParameters().items()}
            if not all(np.isfinite(value) for value in parameters.values()):
                return
            timestep = float(simulation.integrator.getStepSize().value_in_unit(unit.femtosecond))
            self.last = {"xyz": xyz.copy(), "box": box.copy(),
                "step": int(simulation.currentStep), "time_ps": time,
                "potential_kj_mol": energy, "stage": self.stage,
                "context_parameters": parameters, "timestep_fs": timestep,
                "ranking": [{"index": int(i), "force_kj_mol_nm": float(norms[i]),
                    "vector_kj_mol_nm": forces[i].tolist()} for i in indices]}
        except Exception:
            self.problem = "The requested force observation could not be read."

    def failed(self, simulation: Any, error: BaseException) -> None:
        if not self.publish or not numeric_failure(error):
            return
        try:
            self._write(simulation, error)
        except Exception:
            logger.warning("The instability evidence could not be saved; the original failure is unchanged.")
            try:
                target = self.output / "instability.json"
                target.write_text(json.dumps({"version": 1, "available": False,
                    "reason": "A finite force snapshot could not be retained for this failure. No force ranking will be inferred."}), encoding="utf-8")
            except Exception:
                pass

    def _write(self, simulation: Any, error: BaseException) -> None:
        from openmm import app, unit

        last = self.last
        record = {"version": 1, "available": False, "notice": NOTICE,
            "failure": str(error)[:4000], "failure_stage": self.stage,
            "failure_reported_step": int(simulation.currentStep),
            "failure_step_note": "The reported step counter may precede the failed integrator step; the exact failing step is unknown.",
            "interval_steps": self.interval, "reason": self.problem,
            "force_units": "kJ mol^-1 nm^-1", "position_units": "nm",
            "system_sha256": self.model, "integration_force_groups": self.groups}
        record["config_sha256"] = self.config_digest
        if last is not None:
            # This private copy prevents formatting from changing the running
            # topology (PDBFile writes the topology's cell, not the State's).
            from copy import deepcopy

            topology = deepcopy(simulation.topology)
            topology.setPeriodicBoxVectors(last["box"] * unit.nanometer)
            stream = io.StringIO()
            app.PDBFile.writeFile(topology, last["xyz"] * unit.nanometer, stream, keepIds=True)
            snapshot = stream.getvalue()
            raw = snapshot.encode("utf-8")
            lines = [line for line in snapshot.splitlines() if line[:6] in {"ATOM  ", "HETATM"}]
            if len(raw) <= MAX_SNAPSHOT and len(lines) == len(self.atoms):
                digest = hashlib.sha256(raw).hexdigest()
                name = f"instability-snapshot-{digest[:20]}.pdb"
                ranking = []
                for row in last["ranking"]:
                    atom = self.atoms[row["index"]]
                    line = lines[row["index"]]
                    ranking.append({**row,
                        "identity": [line[21:22].strip(), str(int(line[22:26])), line[26:27].strip(),
                                     line[17:20].strip(), line[12:16].strip()],
                        "original_identity": [atom.residue.chain.id, atom.residue.id,
                            atom.residue.insertionCode, atom.residue.name, atom.name],
                        "position_nm": last["xyz"][row["index"]].tolist()})
                (self.output / name).write_bytes(raw)
                record.update({"available": True, "reason": "", "snapshot": name,
                    "snapshot_sha256": digest, "atom_count": len(lines), "ranking": ranking,
                    **{key: value for key, value in last.items() if key not in {"xyz", "box", "ranking"}}})
                if last["stage"] == self.stage and simulation.currentStep >= last["step"]:
                    record["reported_step_gap"] = int(simulation.currentStep) - last["step"]
            else:
                record["reason"] = "The matching full-system snapshot exceeds the display limit."
        target = self.output / "instability.json"
        temporary = target.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(record, indent=2, allow_nan=False), encoding="utf-8")
        temporary.replace(target)
