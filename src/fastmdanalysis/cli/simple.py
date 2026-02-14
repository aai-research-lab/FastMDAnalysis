from __future__ import annotations

from typing import Dict, Any, Callable, Optional
import argparse
import logging

from ._common import add_file_args


def register_simple(subparsers: argparse._SubParsersAction, common_parser: argparse.ArgumentParser) -> None:
    """
    Register simple wrappers around FastMDAnalysis methods. Each entry defines:
    - name: subcommand & method name
    - help: help text
    - args: a function that adds any extra CLI args for that method
    - call: a function that calls the method with parsed args
    """
    specs = [
        {
            "name": "rmsd",
            "help": "RMSD analysis",
            "args": _args_rmsd,
            "call": _call_rmsd,
        },
        {"name": "rmsf", "help": "RMSF analysis", "args": None, "call": _call_passthrough("rmsf")},
        {"name": "rg", "help": "Radius of gyration analysis", "args": None, "call": _call_passthrough("rg")},
        {"name": "hbonds", "help": "Hydrogen bonds analysis", "args": None, "call": _call_passthrough("hbonds")},
        {
            "name": "phi",
            "help": "Phi dihedral analysis",
            "args": _args_dihedral,
            "call": _call_phi,
        },
        {
            "name": "psi",
            "help": "Psi dihedral analysis",
            "args": _args_dihedral,
            "call": _call_psi,
        },
        {
            "name": "omega",
            "help": "Omega dihedral analysis",
            "args": _args_dihedral,
            "call": _call_omega,
        },
        {
            "name": "dihedrals",
            "help": "Combined dihedral analysis (phi/psi/omega + Ramachandran)",
            "args": _args_dihedrals,
            "call": _call_dihedrals,
        },
        {
            "name": "cluster",
            "help": "Clustering analysis",
            "args": _args_cluster,
            "call": _call_cluster,
        },
        {"name": "ss", "help": "Secondary structure (SS) analysis", "args": None, "call": _call_passthrough("ss")},
        {
            "name": "sasa",
            "help": "Solvent accessible surface area (SASA) analysis",
            "args": _args_sasa,
            "call": _call_sasa,
        },
        {
            "name": "dimred",
            "help": "Dimensionality reduction analysis",
            "args": _args_dimred,
            "call": _call_dimred,
        },
        {
            "name": "qvalue",
            "help": "Fraction of native contacts (Q-value) analysis",
            "args": _args_qvalue,
            "call": _call_qvalue,
        },
    ]

    for spec in specs:
        p = subparsers.add_parser(spec["name"], parents=[common_parser], help=spec["help"], conflict_handler="resolve")
        add_file_args(p)
        if spec["args"]:
            spec["args"](p)
        p.set_defaults(_handler=_make_handler(spec["call"], spec["name"]))


def _make_handler(caller: Callable[[Any, argparse.Namespace], Any], name: str):
    def _handler(args: argparse.Namespace, fastmda, logger: logging.Logger) -> None:
        logger.info("Running %s analysis...", name)
        try:
            result = caller(fastmda, args)
            # Try to call .run() if the result exposes it (legacy pattern)
            ran = False
            runner = getattr(result, "run", None)
            if callable(runner):
                result = runner()
                ran = True
            if not ran:
                logger.debug("No .run() method detected; assuming analysis executed inside the method.")
            # Optional plotting
            plotter = getattr(result, "plot", None)
            if callable(plotter):
                plot_res = plotter()
                if isinstance(plot_res, dict):
                    for key, path in plot_res.items():
                        logger.info("Plot for %s saved to: %s", key, path)
                else:
                    logger.info("Plot saved to: %s", plot_res)
            logger.info("%s analysis completed successfully.", name)
        except Exception as e:
            logger.error("Error during %s analysis: %s", name, e)
            raise SystemExit(1)
    return _handler


# --------- Per-method arg adders & callers -----------------------------------

def _args_rmsd(p: argparse.ArgumentParser) -> None:
    # Support: --reference-frame, --ref, and (via argv normalization) -ref
    p.add_argument(
        "--reference-frame", "--ref",
        dest="reference_frame", type=int, default=0,
        help="Reference frame index for RMSD analysis",
    )


def _call_rmsd(fastmda, args: argparse.Namespace):
    return fastmda.rmsd(ref=args.reference_frame, atoms=getattr(args, "atoms", None), output=args.output)


def _args_cluster(p: argparse.ArgumentParser) -> None:
    p.add_argument("--eps", type=float, default=0.5, help="DBSCAN: Maximum distance between samples")
    p.add_argument("--min_samples", type=int, default=5, help="DBSCAN: Minimum samples in a neighborhood")
    p.add_argument("--methods", type=str, nargs="+", default=["dbscan"],
                   help="Clustering methods (e.g., 'dbscan', 'kmeans', 'hierarchical').")
    p.add_argument("--n_clusters", type=int, default=None, help="For KMeans/Hierarchical: number of clusters")


def _call_cluster(fastmda, args: argparse.Namespace):
    return fastmda.cluster(
        methods=args.methods, eps=args.eps, min_samples=args.min_samples,
        n_clusters=args.n_clusters, atoms=getattr(args, "atoms", None), output=args.output
    )


def _args_sasa(p: argparse.ArgumentParser) -> None:
    p.add_argument("--probe_radius", type=float, default=0.14, help="Probe radius (in nm) for SASA calculation")


def _call_sasa(fastmda, args: argparse.Namespace):
    return fastmda.sasa(probe_radius=args.probe_radius, atoms=getattr(args, "atoms", None), output=args.output)


def _args_dimred(p: argparse.ArgumentParser) -> None:
    p.add_argument("--methods", type=str, nargs="+", default=["all"],
                   help="Dimensionality reduction methods (e.g., 'pca', 'mds', 'tsne'). 'all' uses all methods.")


def _call_dimred(fastmda, args: argparse.Namespace):
    return fastmda.dimred(methods=args.methods, atoms=getattr(args, "atoms", None), output=args.output)


def _args_dihedral(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--residues",
        type=int,
        nargs="+",
        default=None,
        help="Residue indices to analyze (0-based). Example: --residues 1 5 10",
    )
    p.add_argument(
        "--units",
        type=str,
        default="degrees",
        choices=["degrees", "radians"],
        help="Units for output angles (default: degrees)",
    )


def _args_dihedrals(p: argparse.ArgumentParser) -> None:
    _args_dihedral(p)
    p.add_argument(
        "--types",
        type=str,
        nargs="+",
        default=["phi", "psi", "omega"],
        help="Dihedral types to compute (phi, psi, omega)",
    )


def _call_phi(fastmda, args: argparse.Namespace):
    return fastmda.phi(
        residues=args.residues,
        units=args.units,
        atoms=getattr(args, "atoms", None),
        output=args.output,
    )


def _call_psi(fastmda, args: argparse.Namespace):
    return fastmda.psi(
        residues=args.residues,
        units=args.units,
        atoms=getattr(args, "atoms", None),
        output=args.output,
    )


def _call_omega(fastmda, args: argparse.Namespace):
    return fastmda.omega(
        residues=args.residues,
        units=args.units,
        atoms=getattr(args, "atoms", None),
        output=args.output,
    )


def _call_dihedrals(fastmda, args: argparse.Namespace):
    return fastmda.dihedrals(
        types=args.types,
        residues=args.residues,
        units=args.units,
        atoms=getattr(args, "atoms", None),
        output=args.output,
    )


def _args_qvalue(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--reference-frame", "--ref",
        dest="reference_frame", type=int, default=0,
        help="Reference frame index for native state (default: 0)",
    )
    p.add_argument(
        "--beta",
        dest="beta_const", type=float, default=50.0,
        help="Beta constant in nm^-1 (default: 50.0)",
    )
    p.add_argument(
        "--lambda",
        dest="lambda_const", type=float, default=1.8,
        help="Lambda constant (default: 1.8)",
    )
    p.add_argument(
        "--cutoff",
        dest="native_cutoff", type=float, default=0.45,
        help="Native contact cutoff distance in nm (default: 0.45)",
    )


def _call_qvalue(fastmda, args: argparse.Namespace):
    return fastmda.qvalue(
        reference_frame=args.reference_frame,
        beta_const=args.beta_const,
        lambda_const=args.lambda_const,
        native_cutoff=args.native_cutoff,
        atoms=getattr(args, "atoms", None),
        output=args.output
    )


def _call_passthrough(method_name: str):
    def _caller(fastmda, args: argparse.Namespace):
        method = getattr(fastmda, method_name)
        return method(atoms=getattr(args, "atoms", None), output=args.output)
    return _caller

