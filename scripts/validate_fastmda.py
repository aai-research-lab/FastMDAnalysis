from __future__ import annotations

import argparse
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import mdtraj as md

from fastmdanalysis import FastMDAnalysis
from fastmdanalysis.analysis.rmsd import _rmsd_no_fit
from fastmdanalysis.analysis.cluster import relabel_compact_positive

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms as mda_rms
    from MDAnalysis.analysis import polymer as mda_polymer
    from MDAnalysis.analysis import sasa as mda_sasa
    from MDAnalysis.analysis import dssp as mda_dssp
except ImportError:
    mda = None
    mda_rms = None
    mda_polymer = None
    mda_sasa = None
    mda_dssp = None

LOG = logging.getLogger("fastmda.validator")


def parse_frames(spec: Optional[str]) -> Optional[Tuple[Optional[int], Optional[int], int]]:

    if spec is None:
        return None
    parts = spec.split(":")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("frames must be formatted as start:stop:step (use blanks for open ends)")

    def _parse(token: str) -> Optional[int]:
        token = token.strip()
        if token == "":
            return None
        return int(token)

    start = _parse(parts[0])
    stop = _parse(parts[1])
    step = _parse(parts[2])
    if step is None or step == 0:
        step = 1
    if step < 0:
        step = -step
    return start, stop, step


def describe_array(arr: np.ndarray) -> Dict[str, Any]:

    arr = np.asarray(arr)
    if arr.size == 0:
        return {"shape": list(arr.shape), "min": None, "max": None, "mean": None, "std": None}
    return {
        "shape": list(arr.shape),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
    }


def compare_arrays(
    fast: np.ndarray,
    reference: np.ndarray,
    tolerance: float,
) -> Dict[str, Any]:

    fast = np.asarray(fast)
    reference = np.asarray(reference)
    info: Dict[str, Any] = {"shape_match": fast.shape == reference.shape}
    if not info["shape_match"]:
        info.update({
            "passed": False,
            "detail": f"shape mismatch fast={fast.shape}, ref={reference.shape}",
        })
        return info
    if fast.size == 0:
        info.update({
            "passed": True,
            "max_abs_diff": 0.0,
            "mean_abs_diff": 0.0,
            "rmse": 0.0,
        })
        return info
    numeric_kinds = {"f", "c", "i", "u", "b"}
    if fast.dtype.kind not in numeric_kinds or reference.dtype.kind not in numeric_kinds:
        mismatches = int(np.sum(fast != reference))
        info.update({
            "passed": mismatches == 0,
            "mismatches": mismatches,
        })
        if mismatches:
            info["detail"] = f"{mismatches} element(s) differ"
        return info
    diff = fast.astype(np.float64) - reference.astype(np.float64)
    info.update({
        "max_abs_diff": float(np.max(np.abs(diff))),
        "mean_abs_diff": float(np.mean(np.abs(diff))),
        "rmse": float(np.sqrt(np.mean(diff ** 2))),
    })
    info["passed"] = info["max_abs_diff"] <= tolerance
    return info


def skipped(reason: str) -> Dict[str, Any]:
    return {"available": False, "detail": reason}


def errored(exc: Exception) -> Dict[str, Any]:
    return {"available": True, "passed": False, "detail": str(exc)}


@dataclass
class MDAnalysisHelper:


    available: bool
    universe: Optional[Any]
    start: Optional[int]
    stop: Optional[int]
    step: int
    reason: Optional[str] = None

    @classmethod
    def build(
        cls,
        top: Path,
        traj: Path,
        frames: Optional[Tuple[Optional[int], Optional[int], int]],
        skip: bool,
    ) -> "MDAnalysisHelper":
        if skip:
            return cls(False, None, None, None, 1, reason="MDAnalysis validation disabled by user")
        if mda is None:
            return cls(False, None, None, None, 1, reason="MDAnalysis not installed")
        try:
            universe = mda.Universe(str(top), str(traj))
        except Exception as exc:
            return cls(False, None, None, None, 1, reason=f"MDAnalysis Universe failed to load: {exc}")
        start = frames[0] if frames else None
        stop = frames[1] if frames else None
        step = frames[2] if frames else 1
        return cls(True, universe, start, stop, step)

    def select(self, selection: Optional[str]):
        if not self.available or self.universe is None:
            raise RuntimeError(self.reason or "MDAnalysis unavailable")
        if selection:
            group = self.universe.select_atoms(selection)
            if group is None or len(group) == 0:
                raise ValueError(f"MDAnalysis selection yielded no atoms: {selection}")
            return group
        return self.universe.atoms

    def frame_kwargs(self) -> Dict[str, Any]:
        return {"start": self.start, "stop": self.stop, "step": self.step}

    def rmsd(self, selection: Optional[str], ref_frame: int) -> np.ndarray:
        group = self.select(selection)
        ref_frame = ref_frame % self.universe.trajectory.n_frames
        analysis = mda_rms.RMSD(
            self.universe,
            self.universe,
            select=selection or "all",
            ref_frame=ref_frame,
            **self.frame_kwargs(),
        )
        analysis.run()
        data = getattr(analysis, "rmsd", getattr(analysis, "results", None))
        if hasattr(data, "rmsd"):
            values = np.asarray(data.rmsd)
        else:
            values = np.asarray(data)
        return values[:, -1] / 10.0

    def rmsf(self, selection: Optional[str]) -> np.ndarray:
        group = self.select(selection)
        analysis = mda_rms.RMSF(group, **self.frame_kwargs())
        analysis.run()
        data = getattr(analysis, "results", None)
        values = np.asarray(getattr(data, "rmsf", getattr(analysis, "rmsf")))
        return values / 10.0

    def radius_of_gyration(self, selection: Optional[str]) -> np.ndarray:
        group = self.select(selection)
        analysis = mda_polymer.RadiusOfGyration(group, **self.frame_kwargs())
        analysis.run()
        data = getattr(analysis, "results", None)
        values = np.asarray(getattr(data, "rg", getattr(analysis, "rg")))
        return values / 10.0

    def sasa(self, selection: Optional[str], probe_radius_nm: float) -> Dict[str, np.ndarray]:
        group = self.select(selection)
        analysis = mda_sasa.ShrakeRupley(
            group,
            probe_radius=probe_radius_nm * 10.0,
            mode="atom",
            **self.frame_kwargs(),
        )
        analysis.run()
        atom_areas = np.asarray(getattr(analysis.results, "sasa"))
        residue_ids = np.array([atom.residue.resindex for atom in group.atoms], dtype=int)
        unique_res = np.unique(residue_ids)
        residue_map = {res: i for i, res in enumerate(unique_res)}
        residue_sasa = np.zeros((atom_areas.shape[0], unique_res.size), dtype=float)
        for res, idx in residue_map.items():
            residue_sasa[:, idx] = atom_areas[:, residue_ids == res].sum(axis=1)
        total = atom_areas.sum(axis=1)
        average = residue_sasa.mean(axis=0)
        scale = 0.01
        return {
            "total": total * scale,
            "residue": residue_sasa * scale,
            "average_per_residue": average * scale,
        }

    def dssp(self, selection: Optional[str]) -> np.ndarray:
        group = self.select(selection)
        analysis = mda_dssp.DSSP(
            self.universe,
            select=selection or "all",
            dssp="mdtraj",
            **self.frame_kwargs(),
        )
        analysis.run()
        data = getattr(analysis, "results", None)
        letters = np.asarray(getattr(data, "secondary_structure", getattr(analysis, "secondary_structure")))
        return letters


@dataclass
class ValidationContext:
    fastmda: FastMDAnalysis
    traj: md.Trajectory
    atoms: Optional[str]
    tolerance: float
    output_root: Path
    reference_frame: int
    align_rmsd: bool
    frames: Optional[Tuple[Optional[int], Optional[int], int]]
    mdanalysis: MDAnalysisHelper

    def atom_indices(self) -> Optional[np.ndarray]:
        if self.atoms is None:
            return None
        sel = self.traj.topology.select(self.atoms)
        if sel is None or len(sel) == 0:
            raise ValueError(f"MDTraj selection yielded no atoms: {self.atoms}")
        return sel

    def subset_traj(self) -> md.Trajectory:
        if self.atoms is None:
            return self.traj
        idx = self.atom_indices()
        return self.traj.atom_slice(idx)


def compute_mdtraj_rmsd(ctx: ValidationContext) -> np.ndarray:
    ref_idx = ctx.reference_frame
    atom_idx = ctx.atom_indices()
    if ctx.align_rmsd:
        values = md.rmsd(ctx.traj, ctx.traj[ref_idx], atom_indices=atom_idx)
    else:
        ref = ctx.traj[ref_idx]
        values = _rmsd_no_fit(ctx.traj, ref, atom_indices=atom_idx)
    return np.asarray(values, dtype=float)


def compute_mdtraj_rmsf(ctx: ValidationContext) -> np.ndarray:
    subtraj = ctx.subset_traj()
    avg_xyz = np.mean(subtraj.xyz, axis=0, keepdims=True)
    ref = md.Trajectory(avg_xyz, subtraj.topology)
    values = md.rmsf(subtraj, ref)
    return np.asarray(values, dtype=float)


def compute_mdtraj_rg(ctx: ValidationContext) -> np.ndarray:
    subtraj = ctx.subset_traj()
    values = md.compute_rg(subtraj)
    return np.asarray(values, dtype=float)


def compute_mdtraj_hbonds(ctx: ValidationContext) -> np.ndarray:
    traj = ctx.traj
    selection_used = ctx.atoms
    if ctx.atoms:
        idx = ctx.traj.topology.select(ctx.atoms)
        if idx is None or len(idx) == 0:
            raise ValueError(f"MDTraj selection yielded no atoms: {ctx.atoms}")
        traj = ctx.traj.atom_slice(idx)
        if traj.topology.n_bonds == 0:
            protein_idx = ctx.traj.topology.select("protein")
            if protein_idx is not None and len(protein_idx) > 0:
                traj = ctx.traj.atom_slice(protein_idx)
                selection_used = "protein"
            else:
                traj = ctx.traj
                selection_used = "all atoms"
    counts = np.zeros(traj.n_frames, dtype=int)
    for i in range(traj.n_frames):
        hb = md.baker_hubbard(traj[i], periodic=False)
        counts[i] = len(hb)
    LOG.debug("HBonds reference used selection: %s", selection_used or "all atoms")
    return counts


def compute_mdtraj_ss(ctx: ValidationContext) -> np.ndarray:
    traj = ctx.subset_traj()
    letters = md.compute_dssp(traj)
    return np.asarray(letters)


def compute_mdtraj_sasa(ctx: ValidationContext, probe_radius: float) -> Dict[str, np.ndarray]:
    traj = ctx.subset_traj()
    try:
        residue_sasa = md.shrake_rupley(traj, probe_radius=probe_radius, mode="residue")
    except TypeError:
        atom_sasa = md.shrake_rupley(traj, probe_radius=probe_radius)
        atom_res = np.array([a.residue.index for a in traj.topology.atoms], dtype=int)
        R = int(atom_res.max() + 1) if atom_res.size else 0
        residue_sasa = np.zeros((traj.n_frames, R), dtype=float)
        for r in range(R):
            mask = atom_res == r
            residue_sasa[:, r] = atom_sasa[:, mask].sum(axis=1)
    total = residue_sasa.sum(axis=1)
    average = residue_sasa.mean(axis=0)
    return {
        "total": np.asarray(total, dtype=float),
        "residue": np.asarray(residue_sasa, dtype=float),
        "average_per_residue": np.asarray(average, dtype=float),
    }


def validate_rmsd(ctx: ValidationContext) -> Dict[str, Any]:
    result: Dict[str, Any] = {"name": "rmsd"}
    analysis = ctx.fastmda.rmsd(
        reference_frame=ctx.reference_frame,
        atoms=ctx.atoms,
        align=ctx.align_rmsd,
        output=str(ctx.output_root / "rmsd_output"),
    )
    fast = np.asarray(analysis.results.get("rmsd"), dtype=float).reshape(-1)
    result["fastmda"] = describe_array(fast)
    comparisons: Dict[str, Any] = {}
    comparisons["mdtraj"] = compare_arrays(fast, compute_mdtraj_rmsd(ctx), ctx.tolerance)
    if not ctx.align_rmsd:
        comparisons["mdanalysis"] = skipped("MDAnalysis RMSD comparison requires alignment (align=True)")
    elif not ctx.mdanalysis.available:
        comparisons["mdanalysis"] = skipped(ctx.mdanalysis.reason or "MDAnalysis unavailable")
    else:
        try:
            comparisons["mdanalysis"] = compare_arrays(fast, ctx.mdanalysis.rmsd(ctx.atoms, ctx.reference_frame), ctx.tolerance)
        except Exception as exc:
            comparisons["mdanalysis"] = errored(exc)
    result["comparisons"] = comparisons
    result["status"] = overall_status(comparisons)
    return result


def validate_rmsf(ctx: ValidationContext) -> Dict[str, Any]:
    result: Dict[str, Any] = {"name": "rmsf"}
    analysis = ctx.fastmda.rmsf(atoms=ctx.atoms, output=str(ctx.output_root / "rmsf_output"))
    fast = np.asarray(analysis.results.get("rmsf"), dtype=float).reshape(-1)
    result["fastmda"] = describe_array(fast)
    comparisons: Dict[str, Any] = {}
    comparisons["mdtraj"] = compare_arrays(fast, compute_mdtraj_rmsf(ctx), ctx.tolerance)
    if not ctx.mdanalysis.available:
        comparisons["mdanalysis"] = skipped(ctx.mdanalysis.reason or "MDAnalysis unavailable")
    else:
        try:
            comparisons["mdanalysis"] = compare_arrays(fast, ctx.mdanalysis.rmsf(ctx.atoms), ctx.tolerance)
        except Exception as exc:
            comparisons["mdanalysis"] = errored(exc)
    result["comparisons"] = comparisons
    result["status"] = overall_status(comparisons)
    return result


def validate_rg(ctx: ValidationContext) -> Dict[str, Any]:
    result: Dict[str, Any] = {"name": "rg"}
    analysis = ctx.fastmda.rg(atoms=ctx.atoms, output=str(ctx.output_root / "rg_output"))
    fast = np.asarray(analysis.results.get("rg"), dtype=float).reshape(-1)
    result["fastmda"] = describe_array(fast)
    comparisons: Dict[str, Any] = {}
    comparisons["mdtraj"] = compare_arrays(fast, compute_mdtraj_rg(ctx), ctx.tolerance)
    if not ctx.mdanalysis.available:
        comparisons["mdanalysis"] = skipped(ctx.mdanalysis.reason or "MDAnalysis unavailable")
    else:
        try:
            comparisons["mdanalysis"] = compare_arrays(fast, ctx.mdanalysis.radius_of_gyration(ctx.atoms), ctx.tolerance)
        except Exception as exc:
            comparisons["mdanalysis"] = errored(exc)
    result["comparisons"] = comparisons
    result["status"] = overall_status(comparisons)
    return result


def validate_hbonds(ctx: ValidationContext) -> Dict[str, Any]:
    result: Dict[str, Any] = {"name": "hbonds"}
    analysis = ctx.fastmda.hbonds(atoms=ctx.atoms, output=str(ctx.output_root / "hbonds_output"))
    fast = np.asarray(analysis.results.get("hbonds_counts"), dtype=float).reshape(-1)
    result["fastmda"] = describe_array(fast)
    comparisons: Dict[str, Any] = {}
    comparisons["mdtraj"] = compare_arrays(fast, compute_mdtraj_hbonds(ctx), ctx.tolerance)
    comparisons["mdanalysis"] = skipped("Dedicated MDAnalysis Baker-Hubbard backend not implemented")
    result["comparisons"] = comparisons
    result["status"] = overall_status(comparisons)
    return result


def validate_ss(ctx: ValidationContext) -> Dict[str, Any]:
    result: Dict[str, Any] = {"name": "ss"}
    analysis = ctx.fastmda.ss(atoms=ctx.atoms, output=str(ctx.output_root / "ss_output"))
    fast_letters = np.asarray(analysis.results.get("ss_data"))
    result["fastmda"] = {"shape": list(fast_letters.shape)}
    comparisons: Dict[str, Any] = {}
    comparisons["mdtraj"] = compare_arrays(fast_letters, compute_mdtraj_ss(ctx), ctx.tolerance)
    if not ctx.mdanalysis.available:
        comparisons["mdanalysis"] = skipped(ctx.mdanalysis.reason or "MDAnalysis unavailable")
    else:
        try:
            comparisons["mdanalysis"] = compare_arrays(fast_letters, ctx.mdanalysis.dssp(ctx.atoms), ctx.tolerance)
        except Exception as exc:
            comparisons["mdanalysis"] = errored(exc)
    result["comparisons"] = comparisons
    result["status"] = overall_status(comparisons)
    return result


def validate_sasa(ctx: ValidationContext, probe_radius: float) -> Dict[str, Any]:
    result: Dict[str, Any] = {"name": "sasa"}
    analysis = ctx.fastmda.sasa(probe_radius=probe_radius, atoms=ctx.atoms, output=str(ctx.output_root / "sasa_output"))
    from_dict: Dict[str, Any] = {}
    if isinstance(analysis.results, dict):
        from_dict.update(analysis.results)
    if isinstance(getattr(analysis, "data", None), dict):
        from_dict.update(analysis.data)
    def _get_array(key: str) -> np.ndarray:
        if key not in from_dict or from_dict.get(key) is None:
            raise KeyError(f"FastMDAnalysis SASA results missing '{key}' data")
        return np.asarray(from_dict[key])
    fast_data = {
        "total": _get_array("total_sasa"),
        "residue": _get_array("residue_sasa"),
        "average_per_residue": _get_array("average_residue_sasa"),
    }
    result["fastmda"] = {k: describe_array(v) for k, v in fast_data.items()}
    ref_mdtraj = compute_mdtraj_sasa(ctx, probe_radius)
    comparisons: Dict[str, Any] = {
        "mdtraj": {
            key: compare_arrays(fast_data[key], ref_mdtraj[key], ctx.tolerance) for key in fast_data
        }
    }
    if not ctx.mdanalysis.available:
        comparisons["mdanalysis"] = {key: skipped(ctx.mdanalysis.reason or "MDAnalysis unavailable") for key in fast_data}
    else:
        try:
            ref_mda = ctx.mdanalysis.sasa(ctx.atoms, probe_radius)
            comparisons["mdanalysis"] = {
                key: compare_arrays(fast_data[key], ref_mda[key], ctx.tolerance) for key in fast_data
            }
        except Exception as exc:
            comparisons["mdanalysis"] = {key: errored(exc) for key in fast_data}
    result["comparisons"] = comparisons
    result["status"] = overall_status_nested(comparisons)
    return result


def flatten_xyz(traj: md.Trajectory, atom_indices: Optional[np.ndarray]) -> np.ndarray:
    if atom_indices is not None:
        data = traj.xyz[:, atom_indices, :]
    else:
        data = traj.xyz
    return data.reshape(traj.n_frames, -1).astype(np.float32)


def validate_dimred(ctx: ValidationContext, methods: Sequence[str]) -> Dict[str, Any]:
    from sklearn.decomposition import PCA
    from sklearn.manifold import MDS, TSNE

    result: Dict[str, Any] = {"name": "dimred"}
    analysis = ctx.fastmda.dimred(methods=methods, atoms=ctx.atoms, outdir=str(ctx.output_root / "dimred_output"))
    fast_results = {k: np.asarray(v) for k, v in analysis.results.items()}
    result["fastmda"] = {k: describe_array(v) for k, v in fast_results.items()}

    atom_idx = ctx.atom_indices()
    X = flatten_xyz(ctx.traj, atom_idx)
    comparisons: Dict[str, Any] = {}
    for method in methods:
        method = method.lower()
        if method == "pca":
            ref = PCA(n_components=2, random_state=analysis.random_state).fit_transform(X)
        elif method == "mds":
            ref = MDS(n_components=2, n_init=4, random_state=analysis.random_state, normalized_stress="auto").fit_transform(X)
        elif method == "tsne":
            perplexity = analysis.tsne_perplexity or max(5, min(30, X.shape[0] // 10 if X.shape[0] > 0 else 5))
            perplexity = min(perplexity, max(1, X.shape[0] - 1))
            ref = TSNE(
                n_components=2,
                perplexity=perplexity,
                max_iter=analysis.tsne_max_iter,
                random_state=analysis.random_state,
                init="pca",
                learning_rate="auto",
            ).fit_transform(X)
        else:
            continue
        comparisons[method] = compare_arrays(fast_results[method], ref.astype(np.float32), ctx.tolerance)
    result["comparisons"] = {"mdtraj": comparisons, "mdanalysis": {m: skipped("Not applicable") for m in comparisons}}
    result["status"] = overall_status(comparisons)
    return result


def compute_dbscan_labels(traj: md.Trajectory, atom_indices: Optional[np.ndarray], eps: float, min_samples: int) -> np.ndarray:
    from sklearn.cluster import DBSCAN

    T = traj.n_frames
    dist = np.empty((T, T), dtype=np.float32)
    for i in range(T):
        if atom_indices is not None:
            dist[:, i] = md.rmsd(traj, traj[i], atom_indices=atom_indices)
        else:
            dist[:, i] = md.rmsd(traj, traj[i])
    dist = 0.5 * (dist + dist.T)
    np.fill_diagonal(dist, 0.0)
    labels_raw = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed").fit_predict(dist)
    labels, _, _ = relabel_compact_positive(labels_raw, start=1, noise_as_last=True)
    return labels


def compute_kmeans_labels(traj: md.Trajectory, atom_indices: Optional[np.ndarray], n_clusters: int, random_state: int) -> np.ndarray:
    from sklearn.cluster import KMeans

    X = flatten_xyz(traj, atom_indices)
    labels = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10).fit_predict(X)
    return labels + 1


def compute_hierarchical_labels(traj: md.Trajectory, atom_indices: Optional[np.ndarray], n_clusters: int) -> np.ndarray:
    from scipy.cluster.hierarchy import linkage, fcluster

    X = flatten_xyz(traj, atom_indices)
    Z = linkage(X, method="ward")
    labels = fcluster(Z, t=n_clusters, criterion="maxclust")
    return labels


def validate_cluster(ctx: ValidationContext, methods: Sequence[str], eps: float, min_samples: int, n_clusters: Optional[int]) -> Dict[str, Any]:
    result: Dict[str, Any] = {"name": "cluster"}
    analysis = ctx.fastmda.cluster(
        methods=methods,
        eps=eps,
        min_samples=min_samples,
        n_clusters=n_clusters,
        atoms=ctx.atoms,
        output=str(ctx.output_root / "cluster_output"),
    )
    fast = analysis.results
    summary = {}
    for method, data in fast.items():
        if isinstance(data, dict) and "labels" in data:
            summary[method] = describe_array(np.asarray(data["labels"]))
    result["fastmda"] = summary

    atom_idx = getattr(analysis, "atom_indices", None)
    if atom_idx is None:
        atom_idx = ctx.atom_indices()
    comparisons: Dict[str, Any] = {}
    for method in methods:
        key = method.lower()
        fast_data = fast.get(key)
        if not isinstance(fast_data, dict) or "labels" not in fast_data:
            continue
        try:
            if key == "dbscan":
                ref_labels = compute_dbscan_labels(ctx.traj, atom_idx, eps=analysis.eps, min_samples=analysis.min_samples)
            elif key == "kmeans":
                ref_labels = compute_kmeans_labels(ctx.traj, atom_idx, n_clusters=analysis.n_clusters or 3, random_state=42)
            elif key == "hierarchical":
                ref_labels = compute_hierarchical_labels(ctx.traj, atom_idx, n_clusters=analysis.n_clusters or 3)
            else:
                continue
            comparisons[key] = compare_arrays(np.asarray(fast_data["labels"]), ref_labels, ctx.tolerance)
        except Exception as exc:
            comparisons[key] = errored(exc)
    result["comparisons"] = {"mdtraj": comparisons, "mdanalysis": {m: skipped("Not available") for m in comparisons}}
    result["status"] = overall_status(comparisons)
    return result


def overall_status(comparisons: Dict[str, Any]) -> str:
    statuses = []
    for info in comparisons.values():
        if isinstance(info, dict) and "passed" in info:
            statuses.append(info["passed"])
    if not statuses:
        return "skipped"
    return "pass" if all(statuses) else "fail"


def overall_status_nested(comparisons: Dict[str, Dict[str, Any]]) -> str:
    statuses: List[bool] = []
    for category in comparisons.values():
        for info in category.values():
            if isinstance(info, dict) and "passed" in info:
                statuses.append(info["passed"])
    if not statuses:
        return "skipped"
    return "pass" if all(statuses) else "fail"


def entry_status(entry: Any) -> str:
    if isinstance(entry, dict):
        if entry.get("available") is False and "passed" not in entry:
            return "skipped"
        if "passed" in entry:
            return "pass" if entry.get("passed") else "fail"
        nested = [entry_status(v) for v in entry.values() if isinstance(v, dict)]
        nested = [s for s in nested if s != "skipped"]
        if not nested:
            return "skipped"
        if all(s == "pass" for s in nested):
            return "pass"
        if any(s == "fail" for s in nested):
            return "fail"
        return "skipped"
    return "skipped"


def entry_summary(entry: Any) -> Tuple[str, str]:
    status = entry_status(entry).upper()
    if not isinstance(entry, dict):
        return status, ""

    if "max_abs_diff" in entry:
        detail = "\n".join(
            [
                f"max={entry['max_abs_diff']:.3e}",
                f"mean={entry.get('mean_abs_diff', 0.0):.3e}",
                f"rmse={entry.get('rmse', 0.0):.3e}",
            ]
        )
        return status, detail

    if "mismatches" in entry:
        return status, f"mismatches={int(entry['mismatches'])}"

    if entry.get("available") is False and entry.get("detail"):
        return status, str(entry["detail"])

    if entry.get("detail") and "passed" in entry:
        return status, str(entry["detail"])

    sub_sections: List[str] = []
    for key, value in entry.items():
        if isinstance(value, dict):
            nested_status, nested_detail = entry_summary(value)
            if not nested_status and not nested_detail:
                continue
            fragment = nested_status
            if nested_detail:
                fragment += " (" + nested_detail.replace("\n", ", ") + ")"
            sub_sections.append(f"{key}: {fragment}")
    detail = "\n".join(sub_sections)
    return status, detail


def render_summary_table(results: Sequence[Dict[str, Any]], output_path: Path, overall: str) -> Optional[Path]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        LOG.warning("Unable to render summary figure: %s", exc)
        return None

    columns = ["Analysis", "MDTraj", "MDTraj detail", "MDAnalysis", "MDAnalysis detail", "Status"]
    rows: List[List[str]] = []
    max_detail_lines = 1
    for res in results:
        comps = res.get("comparisons", {}) if isinstance(res, dict) else {}
        md_status, md_detail = entry_summary(comps.get("mdtraj"))
        mda_status, mda_detail = entry_summary(comps.get("mdanalysis"))
        md_detail = md_detail or "-"
        mda_detail = mda_detail or "-"
        max_detail_lines = max(max_detail_lines, md_detail.count("\n") + 1, mda_detail.count("\n") + 1)
        rows.append([
            str(res.get("name", "?")),
            md_status or entry_status(comps.get("mdtraj")).upper(),
            md_detail,
            mda_status or entry_status(comps.get("mdanalysis")).upper(),
            mda_detail,
            str(res.get("status", "skipped")).upper(),
        ])
    rows.append(["OVERALL", "-", "-", "-", "-", overall.upper()])

    height = max(3.0, 0.6 * len(rows) * max_detail_lines)
    fig, ax = plt.subplots(figsize=(11, height))
    ax.axis("off")
    table = ax.table(cellText=rows, colLabels=columns, loc="center", cellLoc="center")
    try:
        table.auto_set_column_width(col=list(range(len(columns))))
    except Exception:
        pass
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.4)
    fig.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", dpi=200)
    plt.close(fig)
    LOG.info("Wrote summary figure to %s", output_path)
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate FastMDAnalysis outputs against reference implementations.")
    parser.add_argument("--traj", type=Path, default=Path("src/fastmdanalysis/data/trp_cage.dcd"), help="Trajectory file path (e.g., DCD)")
    parser.add_argument("--top", type=Path, default=Path("src/fastmdanalysis/data/trp_cage.pdb"), help="Topology file path (e.g., PDB)")
    parser.add_argument("--frames", type=parse_frames, default=None, help="Frame slice start:stop:step (leave blank fields for open ends)")
    parser.add_argument("--atoms", type=str, default=None, help="MDTraj/MDAnalysis atom selection string (e.g., 'protein and name CA')")
    parser.add_argument("--reference-frame", type=int, default=0, help="Reference frame index for RMSD")
    parser.add_argument("--no-align", action="store_true", help="Disable alignment for RMSD comparisons")
    parser.add_argument("--probe-radius", type=float, default=0.14, help="Probe radius in nm for SASA")
    parser.add_argument("--cluster-eps", type=float, default=0.5, help="DBSCAN epsilon in nm")
    parser.add_argument("--cluster-min-samples", type=int, default=5, help="DBSCAN minimum samples")
    parser.add_argument("--cluster-n-clusters", type=int, default=None, help="Target cluster count for KMeans/Hierarchical (optional)")
    parser.add_argument("--dimred-methods", type=str, default="all", help="Dimensionality reduction methods (comma list or 'all')")
    parser.add_argument("--cluster-methods", type=str, default="all", help="Clustering methods (comma list or 'all')")
    parser.add_argument("--output-root", type=Path, default=Path("validator_outputs"), help="Directory where FastMDAnalysis outputs will be written")
    parser.add_argument("--tolerance", type=float, default=1e-5, help="Absolute tolerance for comparisons")
    parser.add_argument("--skip-mdanalysis", action="store_true", help="Skip MDAnalysis comparisons even if MDAnalysis is installed")
    parser.add_argument("--json", type=Path, default=None, help="Optional path to write JSON summary")
    parser.add_argument("--summary-figure", type=Path, default=None, help="Path to save summary PNG table (default: <output_root>/validation_summary.png)")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    return parser


def resolve_methods(spec: str, default_order: Sequence[str]) -> List[str]:
    if spec.lower() == "all":
        return list(default_order)
    items = [s.strip().lower() for s in spec.split(",") if s.strip()]
    valid = [m for m in default_order if m in items]
    if not valid:
        raise ValueError(f"No valid methods found in '{spec}'. Valid options: {', '.join(default_order)}")
    return valid


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s %(message)s")

    if not args.traj.exists():
        parser.error(f"Trajectory file not found: {args.traj}")
    if not args.top.exists():
        parser.error(f"Topology file not found: {args.top}")

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    fastmda = FastMDAnalysis(str(args.traj), str(args.top), frames=args.frames, atoms=args.atoms)
    ctx = ValidationContext(
        fastmda=fastmda,
        traj=fastmda.traj,
        atoms=args.atoms,
        tolerance=args.tolerance,
        output_root=output_root,
        reference_frame=args.reference_frame,
        align_rmsd=not args.no_align,
        frames=args.frames,
        mdanalysis=MDAnalysisHelper.build(args.top, args.traj, args.frames, args.skip_mdanalysis),
    )

    LOG.info("Running FastMDAnalysis validations")
    results: List[Dict[str, Any]] = []
    validators: List[Callable[[ValidationContext], Dict[str, Any]]] = [
        validate_rmsd,
        validate_rmsf,
        validate_rg,
        validate_hbonds,
        validate_ss,
    ]
    for validator in validators:
        LOG.debug("Validating %s", validator.__name__)
        results.append(validator(ctx))

    results.append(validate_sasa(ctx, args.probe_radius))

    dimred_methods = resolve_methods(args.dimred_methods, ["pca", "mds", "tsne"])
    results.append(validate_dimred(ctx, dimred_methods))

    cluster_methods = resolve_methods(args.cluster_methods, ["dbscan", "kmeans", "hierarchical"])
    results.append(
        validate_cluster(
            ctx,
            cluster_methods,
            eps=args.cluster_eps,
            min_samples=args.cluster_min_samples,
            n_clusters=args.cluster_n_clusters,
        )
    )

    overall = "pass"
    for res in results:
        status = res.get("status", "skipped")
        name = res.get("name")
        LOG.info("%s: %s", name, status)
        comps = res.get("comparisons", {}) if isinstance(res, dict) else {}
        md_status, md_detail = entry_summary(comps.get("mdtraj")) if comps else ("", "")
        mda_status, mda_detail = entry_summary(comps.get("mdanalysis")) if comps else ("", "")
        if md_status or md_detail:
            detail_str = md_detail.replace("\n", " | ") if md_detail else ""
            LOG.info("  MDTraj  -> %s%s", md_status, f" | {detail_str}" if detail_str else "")
        if mda_status or mda_detail:
            detail_str = mda_detail.replace("\n", " | ") if mda_detail else ""
            LOG.info("  MDAnalysis -> %s%s", mda_status, f" | {detail_str}" if detail_str else "")
        if status == "fail":
            overall = "fail"
    LOG.info("Overall validation status: %s", overall)

    summary_path = args.summary_figure or (output_root / "validation_summary.png")
    try:
        render_summary_table(results, summary_path.resolve(), overall)
    except Exception as exc:
        LOG.warning("Failed to write summary figure: %s", exc)

    if args.json:
        payload = {
            "trajectory": str(args.traj),
            "topology": str(args.top),
            "atoms": args.atoms,
            "frames": args.frames,
            "tolerance": args.tolerance,
            "overall": overall,
            "results": results,
        }
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2))
        LOG.info("Wrote JSON report to %s", args.json)

    return 0 if overall == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
