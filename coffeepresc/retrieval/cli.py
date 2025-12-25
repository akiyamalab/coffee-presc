import argparse
import tomllib
from pathlib import Path
from typing import Any

from coffeepresc.models import CoffeePrescSettings as CPSettings
from coffeepresc.retrieval import Database


def add_common_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-v", "--verbosity",
        type=int,
        default=1,
        help="Set verbosity level. (e.g., 0: quiet, 1: normal, 2: verbose)",
    )
    parser.add_argument(
        "--log",
        type=Path,
        help="Path to log file. (will be created if not exists)",
    )
    parser.add_argument(
        "--storage",
        type=Path,
        help="Path to database storage file. (default: data.h5)",
    )
    parser.add_argument(
        "-s", "--setting",
        type=Path,
        help="Path to setting toml file. The above options can be put here.",
    )


def create_command(args: argparse.Namespace) -> None:
    if args.setting is not None:
        with open(args.setting, "rb") as f:
            data = tomllib.load(f)
        params: dict[str, Any] = {
            "verbosity": args.verbosity,
            "log": args.log,
            "storage": args.storage,
            "conformers": args.conformers,
            "fragments": args.fragments,
            "dropped": args.dropped,
            "chunk_size": args.chunk_size,
        }
        for k, v in data.items():
            if k not in params or params[k] is None:
                params[k] = v
        if params.get("log") is None:
            params["log"] = Path("fbdb-create.log")

    db = Database(
        storage_path=Path(params.get("storage")),
        verbosity=params.get("verbosity", 1),
        log_path=Path(params["log"]),
    )
    db.create(
        sdf_path=Path(params["conformers"]),
        fragment_path=Path(params["fragments"]),
        chunk_size=params.get("chunk_size", 100000),
        dropped_path=(
            Path(params["dropped"]) if params.get("dropped") is not None else None
        ),
    )


def search_command(args: argparse.Namespace) -> None:
    cp_settings = CPSettings.parse(path=args.setting, **vars(args))
    cp_settings.validate()

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


def main():
    parser = argparse.ArgumentParser(description="Retrieval System Program")
    subparsers = parser.add_subparsers(dest="command", required=True)

    create_parser = subparsers.add_parser(
        "create", help="Create a new record in the database."
    )
    add_common_argument(create_parser)
    create_parser.add_argument(
        "--conformers",
        type=Path,
        help="The SDF file path to create database (.gz, .bz2 and .xz are also supported)",
    )
    create_parser.add_argument(
        "--fragments",
        type=Path,
        help="Path to representative fragments file. (SDF format)",
    )
    create_parser.add_argument(
        "--dropped",
        type=Path,
        default=None,
        help="The SDF file path to save molecules with errors (default: None)",
    )
    create_parser.add_argument(
        "--chunk_size",
        type=int,
        default=100000,
        help="Maximum Number of molecules per chunk (default: 100000)",
    )
    create_parser.set_defaults(func=create_command)

    search_parser = subparsers.add_parser(
        "search", help="Search compounds in the database."
    )
    add_common_argument(search_parser)
    search_parser.add_argument(
        "--query",
        type=Path,
        help="The query file path. (CSV format)",
    )
    search_parser.add_argument(
        "--matched",
        type=Path,
        help="The output matched compounds file path. (CSV format)",
    )
    search_parser.add_argument(
        "--similarity_th",
        type=float,
        help="Similarity threshold for database search. (default: 0.45)",
    )
    search_parser.set_defaults(func=search_command)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return
    else:
        args.func(args)


if __name__ == "__main__":
    main()
