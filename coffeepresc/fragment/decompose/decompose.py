import os
import subprocess

from coffeepresc.common.typing import StrPath
from coffeepresc.models import DecomposeSettings

__PATH_DECOMPOSE = f"{os.path.dirname(__file__)}/bin/decompose"


def decompose_for_representative(
    ligand_path: StrPath,
    fragment_path: StrPath,
) -> None:
    """
    Decompose molecules into fragments for representative database creation.

    :param ligand_path: Path to the input ligands file (.sdf or .mol2).
    :param fragment_path: Path to the decomposed fragments file (.sdf).
    """
    _decompose_molecules(
        ligand_path=ligand_path,
        fragment_path=fragment_path,
        output_path="/dev/null",
    )

def decompose_for_database(
    ligand_path: StrPath,
    annotated_path: StrPath,
) -> None:
    """
    Decompose molecules into fragments for database creation.

    :param ligand_path: Path to the input ligands file (.sdf or .mol2).
    :param annotated_path: Path to the <fragment_info> annotated ligands file (.sdf).
    """
    _decompose_molecules(
        ligand_path=ligand_path,
        fragment_path="/dev/null",
        output_path=annotated_path,
        ins_fragment_id=True,
    )

def _decompose_molecules(
    ligand_path: StrPath,
    fragment_path: StrPath,
    output_path: StrPath,
    capping_atomic_num: int = -1,
    enable_carbon_capping: bool = False,
    ins_fragment_id: bool = False,
    max_ring_size: int = -1,
    no_merge_solitary: bool = False,
) -> None:
    """
    Decompose molecules into fragments using an external C++ binary.

    :param ligand_path: Path to the input ligands file (.sdf or .mol2).
    :param fragment_path: Path to the decomposed fragments file (.sdf).
    :param output_path: Path to the <fragment_info> annotated ligands file (.sdf).
    :param capping_atomic_num: Atomic number of capping atoms. No capping if set to -1.
    :param enable_carbon_capping: Enable capping even for Carbon atoms.
    :param ins_fragment_id: Enable isotope number injection to mark fragment IDs.
    :param max_ring_size: Maximum ring size. No limit if set to -1.
    :param no_merge_solitary: Disable merging of solitary fragments.
    """
    cmd = [
        __PATH_DECOMPOSE,
        "-l", str(ligand_path),
        "-f", str(fragment_path),
        "-o", str(output_path),
    ]
    cmd += ["--capping_atomic_num", str(capping_atomic_num)]
    if enable_carbon_capping:
        cmd += ["--enable_carbon_capping"]
    if ins_fragment_id:
        cmd += ["--ins_fragment_id"]
    cmd += ["--max_ring_size", str(max_ring_size)]
    if no_merge_solitary:
        cmd += ["--no_merge_solitary"]

    subprocess.run(
        cmd,
        check=True,
    )


class Decompose:
    """
    Orchestator class for fragment decomposition.

    This class provides a high-level interface to decompose molecules into fragments
    """

    def __init__(self, settings=None):
        self.__settings: DecomposeSettings = settings

    @property
    def settings(self) -> DecomposeSettings:
        return self.__settings
    
    def run(
        self,
    ) -> None:
        """
        Run the fragment decomposition process.

        Parameters
        ----------
        ligand_path : StrPath
            Path to the input ligands SDF file.
        annotated_path : StrPath
            Path to the output <fragment_info> annotated ligands SDF file.
        """
        decompose_for_database(
            ligand_path=self.__settings.ligand_path,
            annotated_path=self.__settings.annotated_path,
        )
        print(f"Wrote annotated ligands to {self.__settings.annotated_path}.")
