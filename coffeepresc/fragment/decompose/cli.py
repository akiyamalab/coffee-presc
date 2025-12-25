import argparse
from pathlib import Path

from coffeepresc.fragment.decompose.decompose import decompose_for_database, Decompose
from coffeepresc.models import DecomposeSettings


def main():
    """
    Entry point for fragment decomposition CLI.
    """

    parser = argparse.ArgumentParser(
        description="COFFEE-PRESC Fragment Decomposition Program for FBDB create command"
    )

    io_group = parser.add_argument_group("Input/Output files (required)")
    io_group.add_argument(
        "--decompose_ligands",
        type=Path,
        help="Path to input ligands file. (SDF format)",
    )
    io_group.add_argument(
        "--decompose_annotated",
        type=Path,
        help="Path to output <fragment_info> annotated ligands file. (SDF format)",
    )

    args = parser.parse_args()

    decompose_settings = DecomposeSettings(
        ligand_path=args.decompose_ligands,
        annotated_path=args.decompose_annotated,
    )
    decompose_settings.validate()
    decompose = Decompose(decompose_settings)

    decompose.run()


if __name__ == "__main__":
    main()
