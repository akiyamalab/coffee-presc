import argparse
from pathlib import Path

from coffeepresc.cmpdeval import CmpdEvaluator
from coffeepresc.fragquery import FragQuery
from coffeepresc.models import CoffeePrescSettings as CPSettings
from coffeepresc.retrieval import Database


def main():
    parser = argparse.ArgumentParser(description="COFFEE-PRESC Main Program")

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
        "--output",
        type=Path,
        help="Path to final output file. (CSV format)",
    )

    inter_group = parser.add_argument_group("Intermediate files (optional)")
    inter_group.add_argument(
        "--query",
        type=Path,
        help="Path to query file. (CSV format)",
    )
    inter_group.add_argument(
        "--matched",
        type=Path,
        help="Path to matched compounds file. (CSV format)",
    )
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

    retrieval_group = parser.add_argument_group(
        "Compound retrieval parameters (optional)"
    )
    retrieval_group.add_argument(
        "--storage",
        type=Path,
        help="Path to database storage file. (default: data.h5)",
    )
    retrieval_group.add_argument(
        "--similarity_th",
        type=float,
        help="Similarity threshold for database search. (default: 0.45)",
    )

    score_group = parser.add_argument_group("Compound scoring parameters (optional)")
    score_group.add_argument(
        "--penalty_coef",
        type=float,
        help="Similarity-based penalty coefficient. (default: 8.0)",
    )

    setting_group = parser.add_argument_group("Setting file (optional)")
    setting_group.add_argument(
        "-s", "--setting",
        type=Path,
        help="Path to setting toml file. The above options can be put here.",
    )

    args = parser.parse_args()

    cp_settings = CPSettings.parse(path=args.setting, **vars(args))
    cp_settings.validate()

    fq_settings = cp_settings.fragquery
    fq = FragQuery(fq_settings)
    fq.run(
        receptor_path=cp_settings.receptor_path,
        fragments_path=cp_settings.fragments_path,
        output_path=cp_settings.query_path,
        grid_folder=cp_settings.grid_folder,
        log_path=cp_settings.log_path,
        verbosity=cp_settings.verbosity,
    )

    r_settings = cp_settings.retrieval
    db = Database(
        r_settings.storage_path,
        verbosity=cp_settings.verbosity,
        log_path=cp_settings.log_path,
    )
    db.search(
        query_path=cp_settings.query_path,
        output_path=cp_settings.matched_path,
        similarity_th=r_settings.similarity_th,
    )

    ce_settings = cp_settings.scoring
    ce = CmpdEvaluator(ce_settings)
    ce.run(
        matched_path=cp_settings.matched_path,
        query_path=cp_settings.query_path,
    )
    ce.save(
        output_path=cp_settings.output_path,
        query_row=True,
    )


if __name__ == "__main__":
    main()
