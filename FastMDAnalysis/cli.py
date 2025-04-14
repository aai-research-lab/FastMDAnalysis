"""
Command-line interface for FastMDAnalysis.
Provides subcommands for each analysis, including:
  - rmsd: RMSD analysis.
  - rmsf: RMSF analysis.
  - rg: Radius of gyration analysis.
  - hbonds: Hydrogen bonds analysis.
  - cluster: Clustering analysis.
  - ss: Secondary structure analysis (ss).
  - sasa: Solvent accessible surface area (SASA) analysis.
  - dimred: Dimensionality reduction analysis.
"""

import argparse
import logging
import sys
from pathlib import Path

# Explicitly import analysis modules (import some on demand in dispatch).
from .analysis.rmsd import RMSDAnalysis
from .analysis.rmsf import RMSFAnalysis
from .analysis.rg import RGAnalysis
from .analysis.hbonds import HBondsAnalysis
from .analysis.cluster import ClusterAnalysis
from .analysis.ss import SSAnalysis
from .analysis.sasa import SASAAnalysis
from .utils import load_trajectory

# Configure logging.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(
        description="FastMDAnalysis: Fast Automated MD Trajectory Analysis Using MDTraj"
    )
    subparsers = parser.add_subparsers(dest="command", help="Analysis type", required=True)

    # Helper function: add common arguments.
    def add_common_arguments(subparser):
        subparser.add_argument("-traj", "--trajectory", required=True, help="Path to trajectory file")
        subparser.add_argument("-top", "--topology", required=True, help="Path to topology file")
        subparser.add_argument("-o", "--output", default=None, help="Output directory name")
    
    # RMSD subcommand.
    parser_rmsd = subparsers.add_parser("rmsd", help="RMSD analysis")
    add_common_arguments(parser_rmsd)
    parser_rmsd.add_argument("--ref", type=int, default=0, help="Reference frame index for RMSD")

    # RMSF subcommand.
    parser_rmsf = subparsers.add_parser("rmsf", help="RMSF analysis")
    add_common_arguments(parser_rmsf)
    parser_rmsf.add_argument("--selection", type=str, default="all", help="Atom selection for RMSF analysis")

    # Radius of Gyration subcommand.
    parser_rg = subparsers.add_parser("rg", help="Radius of gyration analysis")
    add_common_arguments(parser_rg)
    
    # Hydrogen Bonds subcommand.
    parser_hbonds = subparsers.add_parser("hbonds", help="Hydrogen bonds analysis")
    add_common_arguments(parser_hbonds)
    
    # Cluster Analysis subcommand.
    parser_cluster = subparsers.add_parser("cluster", help="Clustering analysis")
    add_common_arguments(parser_cluster)
    parser_cluster.add_argument("--eps", type=float, default=0.5, help="DBSCAN: Maximum distance between samples")
    parser_cluster.add_argument("--min_samples", type=int, default=5, help="DBSCAN: Minimum samples in a neighborhood")
    parser_cluster.add_argument("--methods", type=str, nargs='+', default=["dbscan"],
                                help="Clustering methods (e.g., 'dbscan', 'kmeans'). You can list multiple methods.")
    parser_cluster.add_argument("--n_clusters", type=int, default=None, help="For KMeans: number of clusters")

    # Secondary Structure subcommand (renamed to ss).
    parser_ss = subparsers.add_parser("ss", help="Secondary structure analysis (ss)")
    add_common_arguments(parser_ss)

    # SASA subcommand.
    parser_sasa = subparsers.add_parser("sasa", help="Solvent Accessible Surface Area (SASA) analysis")
    add_common_arguments(parser_sasa)
    parser_sasa.add_argument("--probe_radius", type=float, default=0.14, help="Probe radius (in nm) for SASA calculation")

    # Dimensionality Reduction subcommand.
    parser_dimred = subparsers.add_parser("dimred", help="Dimensionality reduction analysis")
    add_common_arguments(parser_dimred)
    parser_dimred.add_argument("--methods", type=str, nargs='+', default=["all"],
                               help="Dimensionality reduction methods (e.g., 'pca', 'mds', 'tsne'). 'all' uses all methods.")
    parser_dimred.add_argument("--atom_selection", type=str, default="protein and name CA",
                               help="Atom selection for constructing the feature matrix.")
    
    return parser.parse_args()

def main():
    args = parse_args()

    # Load the trajectory.
    try:
        traj = load_trajectory(args.trajectory, args.topology)
        logger.info(f"Loaded trajectory with {traj.n_frames} frames")
    except Exception as e:
        logger.error(f"Failed to load trajectory: {e}")
        sys.exit(1)

    # Prepare output directory (default is {command}_output if not provided).
    output_dir = args.output or f"{args.command}_output"

    try:
        if args.command == "rmsd":
            analysis = RMSDAnalysis(traj, ref_frame=args.ref, output=output_dir)
        elif args.command == "rmsf":
            analysis = RMSFAnalysis(traj, selection=args.selection, output=output_dir)
        elif args.command == "rg":
            analysis = RGAnalysis(traj, output=output_dir)
        elif args.command == "hbonds":
            analysis = HBondsAnalysis(traj, output=output_dir)
        elif args.command == "cluster":
            analysis = ClusterAnalysis(traj, eps=args.eps, min_samples=args.min_samples,
                                       methods=args.methods, n_clusters=args.n_clusters, output=output_dir)
        elif args.command == "ss":
            analysis = SSAnalysis(traj, output=output_dir)
        elif args.command == "sasa":
            analysis = SASAAnalysis(traj, probe_radius=args.probe_radius, output=output_dir)
        elif args.command == "dimred":
            from .analysis.dimred import DimRedAnalysis
            analysis = DimRedAnalysis(traj, methods=args.methods, atom_selection=args.atom_selection, output=output_dir)
        else:
            logger.error("Unknown command")
            sys.exit(1)

        logger.info(f"Running {args.command} analysis...")
        analysis.run()
        logger.info(f"{args.command} analysis completed successfully.")
        if hasattr(analysis, "plot") and callable(analysis.plot):
            plot_path = analysis.plot()
            logger.info(f"Plot saved to {plot_path}")
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

