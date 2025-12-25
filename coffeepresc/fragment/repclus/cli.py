import argparse
from pathlib import Path

from coffeepresc.fragment.repclus.clustering import Clustering
from coffeepresc.models import RepclusSettings


def main():
    """
    Entry point for representative fragment clustering CLI.
    """

    parser = argparse.ArgumentParser(
        description="COFFEE-PRESC Cluster Representative Fragment Generation Program"
    )

    io_group = parser.add_argument_group("Input/Output files (required)")
    io_group.add_argument(
        "--clustering_molecules",
        type=Path,
        help="Path to input molecules file. (SDF format)",
    )
    io_group.add_argument(
        "--clustering_output",
        type=Path,
        help="Path to output representative fragments file. (SDF format)",
    )

    inner_group = parser.add_argument_group("Clustering parameters (required)")
    inner_group.add_argument(
        "--n_clusters",
        type=int,
        help="Number of clusters to generate.",
    )
    
    args = parser.parse_args()

    clustering_settings = RepclusSettings(
        clustering_molecules=args.clustering_molecules,
        clustering_output=args.clustering_output,
        n_clusters=args.n_clusters,
    )
    clustering_settings.validate()
    clustering = Clustering(clustering_settings)

    clustering.run()


if __name__ == "__main__":
    main()
