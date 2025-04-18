#!/usr/bin/env python
"""
Command-line interface for FastMDAnalysis.
Provides subcommands for various MD analyses:
  - rmsd: RMSD analysis.
  - rmsf: RMSF analysis.
  - rg: Radius of gyration analysis.
  - hbonds: Hydrogen bonds analysis.
  - cluster: Clustering analysis.
  - ss: Secondary structure (SS) analysis.
  - sasa: Solvent accessible surface area (SASA) analysis.
  - dimred: Dimensionality reduction analysis.

Global options:
  --frames  : Frame selection as an iterable "start,stop,stride" (e.g., "0,-1,10"). Negative indices are allowed.
  --atoms   : Global atom selection string (e.g., "protein", "protein and name CA").
  --verbose : When specified, print detailed INFO messages to the screen.
  
File-related options (-traj, -top, -o) are provided at the subcommand level.
"""

import argparse
import logging
import sys
from pathlib import Path

# Create a parent parser for global arguments.
common_parser = argparse.ArgumentParser(add_help=False)
common_parser.add_argument("--frames", type=str, default=None,
                           help="Frame selection as 'start,stop,stride' (e.g., '0,-1,10'). Negative indices allowed.")
common_parser.add_argument("--atoms", type=str, default=None,
                           help='Global atom selection string (e.g., "protein", "protein and name CA").')
common_parser.add_argument("--verbose", action="store_true",
                           help="Print detailed INFO messages to the screen.")

# Helper function: add file-related arguments to each subcommand.
def add_file_args(subparser):
    subparser.add_argument("-traj", "--trajectory", required=True, help="Path to trajectory file")
    subparser.add_argument("-top", "--topology", required=True, help="Path to topology file")
    subparser.add_argument("-o", "--output", default=None, help="Output directory name")

# Main parser including global options.
parser = argparse.ArgumentParser(
    description="FastMDAnalysis: Fast Automated MD Trajectory Analysis Using MDTraj",
    parents=[common_parser]
)
subparsers = parser.add_subparsers(dest="command", help="Analysis type", required=True)

# Subcommand: RMSD.
parser_rmsd = subparsers.add_parser("rmsd", parents=[common_parser], help="RMSD analysis")
add_file_args(parser_rmsd)
parser_rmsd.add_argument("--ref", type=int, default=0, help="Reference frame index for RMSD")
parser_rmsd.add_argument("--selection", type=str, default=None,
                         help="Atom selection for RMSD analysis (overrides global --atoms)")

# Subcommand: RMSF.
parser_rmsf = subparsers.add_parser("rmsf", parents=[common_parser], help="RMSF analysis")
add_file_args(parser_rmsf)
parser_rmsf.add_argument("--selection", type=str, default=None,
                         help="Atom selection for RMSF analysis (overrides global --atoms)")

# Subcommand: RG.
parser_rg = subparsers.add_parser("rg", parents=[common_parser], help="Radius of gyration analysis")
add_file_args(parser_rg)

# Subcommand: HBonds.
parser_hbonds = subparsers.add_parser("hbonds", parents=[common_parser], help="Hydrogen bonds analysis")
add_file_args(parser_hbonds)

# Subcommand: Cluster.
parser_cluster = subparsers.add_parser("cluster", parents=[common_parser], help="Clustering analysis")
add_file_args(parser_cluster)
parser_cluster.add_argument("--eps", type=float, default=0.5, help="DBSCAN: Maximum distance between samples")
parser_cluster.add_argument("--min_samples", type=int, default=5, help="DBSCAN: Minimum samples in a neighborhood")
parser_cluster.add_argument("--methods", type=str, nargs='+', default=["dbscan"],
                            help="Clustering methods (e.g., 'dbscan', 'kmeans').")
parser_cluster.add_argument("--n_clusters", type=int, default=None, help="For KMeans: number of clusters")

# Subcommand: SS.
parser_ss = subparsers.add_parser("ss", parents=[common_parser], help="Secondary structure (SS) analysis")
add_file_args(parser_ss)

# Subcommand: SASA.
parser_sasa = subparsers.add_parser("sasa", parents=[common_parser], help="Solvent accessible surface area (SASA) analysis")
add_file_args(parser_sasa)
parser_sasa.add_argument("--probe_radius", type=float, default=0.14, help="Probe radius (in nm) for SASA calculation")

# Subcommand: Dimensionality Reduction.
parser_dimred = subparsers.add_parser("dimred", parents=[common_parser], help="Dimensionality reduction analysis")
add_file_args(parser_dimred)
parser_dimred.add_argument("--methods", type=str, nargs='+', default=["all"],
                           help="Dimensionality reduction methods (e.g., 'pca', 'mds', 'tsne'). 'all' uses all methods.")
parser_dimred.add_argument("--atom_selection", type=str, default=None,
                           help="Atom selection for constructing the feature matrix (overrides global --atoms)")

def main():
    args = parser.parse_args()

    # Set up logging. Create a log file based on the command name in the output directory.
    output_dir = args.output or f"{args.command}_output"
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    log_filename = Path(output_dir) / f"{args.command}.log"

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers = []  # Reset handlers if any exist.

    # File handler: always log INFO-level details.
    fh = logging.FileHandler(log_filename)
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # Stream handler: log detailed INFO if --verbose is set; otherwise, only WARNING or higher.
    sh = logging.StreamHandler()
    if args.verbose:
        sh.setLevel(logging.INFO)
    else:
        sh.setLevel(logging.WARNING)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    logger.info("FastMDAnalysis command: %s", " ".join(sys.argv))
    logger.info("Parsed arguments: %s", args)

    # Parse the frames argument (expecting a comma-separated "start,stop,stride").
    frames = None
    if args.frames:
        try:
            frames = tuple(map(int, args.frames.split(',')))
            if len(frames) != 3:
                raise ValueError
        except ValueError:
            logger.error("Invalid --frames format. Expected 'start,stop,stride' (e.g., '0,-1,10').")
            sys.exit(1)

    # Use global atoms option.
    atoms = args.atoms

    # Initialize FastMDAnalysis.
    from . import FastMDAnalysis
    fastmda = FastMDAnalysis(args.trajectory, args.topology, frames=frames, atoms=atoms)

    try:
        # Dispatch to the appropriate analysis method based on the subcommand.
        if args.command == "rmsd":
            analysis = fastmda.rmsd(ref=args.ref, atoms=args.selection)
        elif args.command == "rmsf":
            analysis = fastmda.rmsf(atoms=args.selection)
        elif args.command == "rg":
            analysis = fastmda.rg()
        elif args.command == "hbonds":
            analysis = fastmda.hbonds()
        elif args.command == "cluster":
            analysis = fastmda.cluster(methods=args.methods, eps=args.eps,
                                       min_samples=args.min_samples, n_clusters=args.n_clusters)
        elif args.command == "ss":
            analysis = fastmda.ss()
        elif args.command == "sasa":
            analysis = fastmda.sasa(probe_radius=args.probe_radius)
        elif args.command == "dimred":
            analysis = fastmda.dimred(methods=args.methods, atom_selection=args.atom_selection)
        else:
            logger.error("Unknown command")
            sys.exit(1)

        logger.info("Running %s analysis...", args.command)
        analysis.run()
        logger.info("%s analysis completed successfully.", args.command)

        # If analysis supports plotting, generate plots.
        if hasattr(analysis, "plot") and callable(analysis.plot):
            plot_result = analysis.plot()
            if isinstance(plot_result, dict):
                for key, path in plot_result.items():
                    logger.info("Plot for %s saved to: %s", key, path)
            else:
                logger.info("Plot saved to: %s", plot_result)
    except Exception as e:
        logger.error("Analysis failed: %s", e)
        sys.exit(1)

if __name__ == "__main__":
    main()

