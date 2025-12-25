import argparse
from pathlib import Path

from coffeepresc.fragquery import FragQuery
from coffeepresc.models import CoffeePrescSettings as CPSettings


def main():
    parser = argparse.ArgumentParser(
        description="COFFEE-PRESC Fragment Docking & Query Enumeration Program"
    )

    io_group = parser.add_argument_group("Input/Output files (required)")
    io_group.add_argument(
        "--docking_config",
        type=Path,
        help="Path to docking region configuration file.",
    )
    io_group.add_argument(
        "--receptor",
        type=Path,
        help="Path to receptor file. (PDB format)",
    )
    io_group.add_argument(
        "--fragments",
        type=Path,
        help="Path to representative fragments file. (SDF format)",
    )
    io_group.add_argument(
        "--query",
        type=Path,
        help="Path to query file. (CSV format)",
    )

    inter_group = parser.add_argument_group("Intermediate files (optional)")
    inter_group.add_argument(
        "--grid",
        type=Path,
        help="Path to grid folder. (will be created if not exists)",
    )

    log_group = parser.add_argument_group("Logging settings (optional)")
    log_group.add_argument(
        "-v", "--verbosity",
        type=int,
        help="Set verbosity level. (e.g., 0: quiet, 1: normal, 2: verbose)",
    )
    log_group.add_argument(
        "--log",
        type=Path,
        help="Path to log file. (will be created if not exists)",
    )

    query_group = parser.add_argument_group("Query enumeration parameters (optional)")
    query_group.add_argument(
        "--candidate_poses",
        type=int,
        help="Number of candidate poses to be stored for each fragment. (default: 40)",
    )
    query_group.add_argument(
        "--cluster_size",
        type=float,
        help="Clustering size for fragment poses. (default: 1.0)",
    )
    query_group.add_argument(
        "--distance_width",
        type=float,
        help="Distance width for scoring. (default: 0.1)",
    )
    query_group.add_argument(
        "--no_oberrorlog",
        action="store_true",
        default=None,
        help="If set, OpenBabel error messages will not be logged.",
    )

    setting_group = parser.add_argument_group("Setting file (optional)")
    setting_group.add_argument(
        "-s", "--setting",
        type=Path,
        help="Path to setting toml file. The above options can be put here.",
    )

    args = parser.parse_args()

    cp_settings = CPSettings.parse(path=args.setting, **vars(args))

    fq_settings = cp_settings.fragquery
    fq_settings.validate()
    fq = FragQuery(fq_settings)
    fq.run(
        receptor_path=cp_settings.receptor_path,
        fragments_path=cp_settings.fragments_path,
        output_path=cp_settings.query_path,
        grid_folder=cp_settings.grid_folder,
        log_path=cp_settings.log_path,
        verbosity=cp_settings.verbosity,
    )


if __name__ == "__main__":
    main()
