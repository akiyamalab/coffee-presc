import argparse
from pathlib import Path

from coffeepresc.cmpdeval import CmpdEvaluator
from coffeepresc.models import CoffeePrescSettings as CPSettings


def main():
    parser = argparse.ArgumentParser(description="COFFEE-PRESC Main Program")

    io_group = parser.add_argument_group("Input/Output files (required)")
    io_group.add_argument(
        "--query",
        type=Path,
        help="Path to query file. (CSV format)",
    )
    io_group.add_argument(
        "--matched",
        type=Path,
        help="Path to matched compounds file. (CSV format)",
    )
    io_group.add_argument(
        "--output",
        type=Path,
        help="Path to final output file. (CSV format)",
    )

    # log_group = parser.add_argument_group("Logging settings (optional)")
    # log_group.add_argument(
    #     "-v", "--verbosity",
    #     type=int,
    #     help="Set verbosity level. (e.g., 0: quiet, 1: normal, 2: verbose)",
    # )
    # log_group.add_argument(
    #     "--log",
    #     type=Path,
    #     help="Path to log file. (will be created if not exists)",
    # )

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
        help="Path to setting toml file. The above options can be put here."
    )

    args = parser.parse_args()

    cp_settings = CPSettings.parse(path=args.setting, **vars(args))

    ce_settings = cp_settings.scoring
    ce_settings.validate()
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
