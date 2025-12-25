from pathlib import Path
import gzip

from coffeepresc.common.typing import StrPath
from rdkit import Chem


def read_sdf(
    path: StrPath,
    sanitize: bool = True,
    remove_none: bool = True,
) -> list[Chem.Mol]:
    """
    Read molecules from an SDF file robustly.

    This function supports both `.sdf` and compressed `.sdf.gz` files.
    Invalid molecules (returned as `None` by RDKit) are optionally removed.

    Parameters
    ----------
    path : StrPath
        Path to the input SDF file.
    sanitize : bool, default=True
        Whether to sanitize molecules during loading.
    remove_none : bool, default=True
        Whether to remove invalid molecules (`None`) from the result.

    Returns
    -------
    list[rdkit.Chem.Mol]
        List of successfully loaded RDKit molecule objects.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    ValueError
        If no valid molecules are found in the file.
    """

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"SDF file not found: {path}")

    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fin:
            suppl = Chem.ForwardSDMolSupplier(fin, sanitize=sanitize, removeHs=False)
    else:  # path.suffix in [".sdf"]:
        suppl = Chem.ForwardSDMolSupplier(str(path), sanitize=sanitize, removeHs=False)

    molecules: list[Chem.Mol] = [
        molecule
        for molecule in suppl
        if (molecule is not None or not remove_none)        
    ]

    if len(molecules) == 0:
        raise ValueError(f"No valid molecules found in SDF file: {path}")

    return molecules


def write_sdf(
    molecules: list[Chem.Mol],
    path: StrPath,
    overwrite: bool = True,
) -> None:
    """
    Write molecules to an SDF file safely.

    The output directory is created automatically if it does not exist.

    Parameters
    ----------
    molecules : list[rdkit.Chem.Mol]
        Molecules to write.
    path : StrPath
        Path to the output SDF file.
    overwrite : bool, default=True
        Whether to overwrite an existing file.

    Raises
    ------
    FileExistsError
        If the output file exists and `overwrite` is False.
    ValueError
        If the molecule list is empty.
    """

    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"SDF file already exists: {path}")

    path.parent.mkdir(parents=True, exist_ok=True)

    molecules = list(molecules)
    if len(molecules) == 0:
        raise ValueError("No molecules to write to SDF file.")

    writer = Chem.SDWriter(str(path))
    try:
        for molecule in molecules:
            writer.write(molecule)
    finally:
        writer.close()
