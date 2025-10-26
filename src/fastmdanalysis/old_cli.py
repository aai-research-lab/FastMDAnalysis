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
  --frames  : Frame selection as a comma-separated "start,stop,stride" (e.g., "0,-1,10"). Negative indices are allowed.
  --atoms   : Global atom selection string (e.g., "protein", "protein and name CA").
  --verbose : When specified, print detailed log messages (DEBUG and INFO) to the screen.
  
File-related options (-traj, -top, -o) are provided at the subcommand level.
"""

import sys
import argparse
import logging
from pathlib import Path

# Parent parser for global options.
common_parser = argparse.ArgumentParser(add_help=False)
common_parser.add_argument("--frames", type=str, default=None,
                           help="Frame selection as 'start,stop,stride' (e.g., '0,-1,10'). Negative indices allowed.")
common_parser.add_argument("--atoms", type=str, default=None,
                           help='Global atom selection string (e.g., "protein", "protein and name CA").')
common_parser.add_argument("--verbose", action="store_true",
                           help="Print detailed log messages to the screen.")

def add_file_args(subparser):
    subparser.add_argument("-traj", "--trajectory", required=True, help="Path to trajectory file")
    subparser.add_argument("-top", "--topology", required=True, help="Path to topology file")
    subparser.add_argument("-o", "--output", default=None, help="Output directory name")

# Main parser.
parser = argparse.ArgumentParser(
    description="FastMDAnalysis: Fast Automated MD Trajectory Analysis",
    epilog="Docs: https://fastmdanalysis.readthedocs.io/en/latest/", 
    parents=[common_parser]
)
subparsers = parser.add_subparsers(dest="command", help="Analysis type", required=True)

# Subcommand: RMSD.
parser_rmsd = subparsers.add_parser("rmsd", parents=[common_parser], help="RMSD analysis", conflict_handler="resolve")
add_file_args(parser_rmsd)
parser_rmsd.add_argument(
    "--reference-frame",
    dest="reference_frame",
    type=int,
    default=0,
    help="Reference frame index for RMSD analysis",
)

# Subcommand: RMSF.
parser_rmsf = subparsers.add_parser("rmsf", parents=[common_parser], help="RMSF analysis", conflict_handler="resolve")
add_file_args(parser_rmsf)

# Subcommand: RG.
parser_rg = subparsers.add_parser("rg", parents=[common_parser], help="Radius of gyration analysis", conflict_handler="resolve")
add_file_args(parser_rg)

# Subcommand: HBonds.
parser_hbonds = subparsers.add_parser("hbonds", parents=[common_parser], help="Hydrogen bonds analysis", conflict_handler="resolve")
add_file_args(parser_hbonds)

# Subcommand: Cluster.
parser_cluster = subparsers.add_parser("cluster", parents=[common_parser], help="Clustering analysis", conflict_handler="resolve")
add_file_args(parser_cluster)
parser_cluster.add_argument("--eps", type=float, default=0.5, help="DBSCAN: Maximum distance between samples")
parser_cluster.add_argument("--min_samples", type=int, default=5, help="DBSCAN: Minimum samples in a neighborhood")
parser_cluster.add_argument("--methods", type=str, nargs='+', default=["dbscan"],
                            help="Clustering methods (e.g., 'dbscan', 'kmeans', 'hierarchical').")
parser_cluster.add_argument("--n_clusters", type=int, default=None, help="For KMeans/Hierarchical: number of clusters")

# Subcommand: SS.
parser_ss = subparsers.add_parser("ss", parents=[common_parser], help="Secondary structure (SS) analysis", conflict_handler="resolve")
add_file_args(parser_ss)

# Subcommand: SASA.
parser_sasa = subparsers.add_parser("sasa", parents=[common_parser], help="Solvent accessible surface area (SASA) analysis", conflict_handler="resolve")
add_file_args(parser_sasa)
parser_sasa.add_argument("--probe_radius", type=float, default=0.14, help="Probe radius (in nm) for SASA calculation")

# Subcommand: Dimensionality Reduction.
parser_dimred = subparsers.add_parser("dimred", parents=[common_parser], help="Dimensionality reduction analysis", conflict_handler="resolve")
add_file_args(parser_dimred)
parser_dimred.add_argument("--methods", type=str, nargs='+', default=["all"],
                           help="Dimensionality reduction methods (e.g., 'pca', 'mds', 'tsne'). 'all' uses all methods.")

def main():
    args = parser.parse_args()

    # Determine the output directory.
    output_dir = args.output if args.output else f"{args.command}_output"
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    log_filename = Path(output_dir) / f"{args.command}.log"

    # Configure logging using basicConfig.
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_filename, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger()
    logger.info("FastMDAnalysis command: %s", " ".join(sys.argv))
    logger.info("Parsed arguments: %s", args)

    # Parse the frames argument.
    frames = None
    if args.frames:
        try:
            frames = tuple(map(int, args.frames.split(',')))
            if len(frames) != 3:
                raise ValueError
        except ValueError:
            logger.error("Invalid --frames format. Expected 'start,stop,stride' (e.g., '0,-1,10').")
            sys.exit(1)

    # Global atom selection.
    atoms = getattr(args, "atoms", None)

    # Initialize FastMDAnalysis instance (do not pass output here).
    try:
        from fastmdanalysis import FastMDAnalysis
        fastmda = FastMDAnalysis(args.trajectory, args.topology, frames=frames, atoms=atoms)
    except Exception as e:
        logger.error("Error initializing FastMDAnalysis: %s", e)
        sys.exit(1)

    # Dispatch to appropriate analysis.
    try:
        if args.command == "rmsd":
            analysis = fastmda.rmsd(reference_frame=args.reference_frame, atoms=atoms, output=args.output)
        elif args.command == "rmsf":
            analysis = fastmda.rmsf(atoms=atoms, output=args.output)
        elif args.command == "rg":
            analysis = fastmda.rg(atoms=atoms, output=args.output)
        elif args.command == "hbonds":
            analysis = fastmda.hbonds(atoms=atoms, output=args.output)
        elif args.command == "cluster":
            analysis = fastmda.cluster(methods=args.methods, eps=args.eps,
                                       min_samples=args.min_samples, n_clusters=args.n_clusters,
                                       atoms=atoms, output=args.output)
        elif args.command == "ss":
            analysis = fastmda.ss(atoms=atoms, output=args.output)
        elif args.command == "sasa":
            analysis = fastmda.sasa(probe_radius=args.probe_radius, atoms=atoms, output=args.output)
        elif args.command == "dimred":
            analysis = fastmda.dimred(methods=args.methods, atoms=atoms, output=args.output)
        else:
            logger.error("Unknown command: %s", args.command)
            sys.exit(1)

        logger.info("Running %s analysis...", args.command)
        analysis.run()
        logger.info("%s analysis completed successfully.", args.command)

        if hasattr(analysis, "plot") and callable(analysis.plot):
            plot_result = analysis.plot()
            if isinstance(plot_result, dict):
                for key, path in plot_result.items():
                    logger.info("Plot for %s saved to: %s", key, path)
            else:
                logger.info("Plot saved to: %s", plot_result)
    except Exception as e:
        logger.error("Error during %s analysis: %s", args.command, e)
        sys.exit(1)

if __name__ == "__main__":
    main()

