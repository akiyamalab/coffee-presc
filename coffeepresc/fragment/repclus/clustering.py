import tempfile

from coffeepresc.common.sdf_io import read_sdf, write_sdf
from coffeepresc.common.typing import StrPath
from coffeepresc.fragment.decompose import decompose_for_representative
from coffeepresc.fragment.repclus.clusters import Clusters
from coffeepresc.fragment.repclus.filtering import Filter
from coffeepresc.models import RepclusSettings
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.warning")
RDLogger.DisableLog("rdApp.error")


def _build_molecule_filter() -> Filter:
    """
    Construct and return a molecule-level filter.

    The filter is configured with predefined ligand-level filtering rules
    suiatable for representative fragment generation.

    Returns
    -------
    Filter
        A configured Filter instance for molecule preprocessing.
    """

    molecule_filter: Filter = Filter()
    molecule_filter.set_ligand_filters()
    return molecule_filter


def _build_fragment_filter() -> Filter:
    """
    Construct and return a fragment-level filter.

    The filter is configured with predefined fragment-level filtering rules
    suitable for representative fragment generation.

    Returns
    -------
    Filter
        A configured Filter instance for fragment preprocessing.
    """

    fragment_filter: Filter = Filter()
    fragment_filter.set_fragment_filters()
    return fragment_filter


def _clustering_flow(
    molecules_path: StrPath,
    output_representative_path: StrPath,
    tmp_molecules_path: StrPath,
    tmp_fragments_path: StrPath,
    molecule_filter: Filter,
    fragment_filter: Filter,
    n_clusters: int,
) -> dict[str, int]:
    """
    Execute the full clustering workflow for representative fragment selection.

    This function performs the following steps:
    1. Load input molecules from the SDF file.
    2. Apply molecule-level filtering.
    3. Decompose filtered molecules into fragments.
    4. Apply fragment-level filtering.
    5. Perform Ward hierarchical clustering on the fragments.
    6. Select representative fragments from each cluster.
    7. Write the representative fragments to the output SDF file.

    Parameters
    ----------
    molecules_path : StrPath
        Path to the input SDF file containing molecules.
    output_representative_path : StrPath
        Path to the output SDF file for representative fragments.
    tmp_molecules_path : StrPath
        Path to a temporary SDF file for filtered molecules.
    tmp_fragments_path : StrPath
        Path to a temporary SDF file for filtered fragments.
    molecule_filter : Filter
        Filter instance for molecule-level preprocessing.
    fragment_filter : Filter
        Filter instance for fragment-level preprocessing.
    n_clusters : int
        Number of clusters to form during clustering.

    Returns
    -------
    dict[str, int]
        Statistics of the clustering process, including:
        - number of input molecules,
        - number of filtered molecules,
        - number of generated fragments,
        - number of filtered fragments,
        - number of representative fragments.

    Raise
    -------
    ValueError
        If no molecules or fragments pass the filtering steps.
    """

    # 1. Load and filter molecules
    print("Loading molecules...")
    molecules: list[Chem.Mol] = read_sdf(molecules_path)
    print(f"Loaded {len(molecules)} molecules.")

    # 2. Apply molecule-level filtering
    filtered_molecules: list[Chem.Mol] = molecule_filter.apply(molecules)
    if not filtered_molecules:
        raise ValueError("No molecules passed the filtering step.")
    print(f"{len(filtered_molecules)} molecules passed filtering.")
    write_sdf(filtered_molecules, tmp_molecules_path, overwrite=True)

    # 3. Decompose molecules into fragments and save to temporary file
    print("Decomposing molecules into fragments...")
    decompose_for_representative(tmp_molecules_path, tmp_fragments_path)
    print("Decomposition completed.")

    # 4. Load and filter fragments
    print("Loading fragments...")
    fragments: list[Chem.Mol] = read_sdf(tmp_fragments_path)
    filtered_fragments: list[Chem.Mol] = fragment_filter.apply(fragments)
    if not filtered_fragments:
        raise ValueError("No fragments passed the filtering step.")
    print(f"{len(filtered_fragments)} fragments passed filtering.")

    # 4.5 tagged fragments
    for i, frag in enumerate(filtered_fragments):
        frag.SetProp("SMILES", Chem.MolToSmiles(frag))
        frag.SetProp("fragment_id", f"{Chem.MolToSmiles(frag)}_{i+1}")

    # 5. Perform clustering
    print("Clustering fragments...")
    clusters: Clusters = Clusters.from_fragments_by_ward(
        molecules=filtered_fragments,
        n_cluster=n_clusters,
    )
    print(f"Formed {len(clusters.clusters)} clusters.")

    # 6. Build representative fragments
    print("Selecting representative fragments...")
    clusters.build_representatives()
    representatives: list[Chem.Mol] = clusters.get_representatives()
    print(f"Selected {len(representatives)} representative fragments.")

    # 7. Write representative fragments to output file
    write_sdf(representatives, output_representative_path, overwrite=True)
    print(f"Wrote representative fragments to {output_representative_path}.")

    return {
        "num_input_molecules": len(molecules),
        "num_filtered_molecules": len(filtered_molecules),
        "num_input_fragments": len(fragments),
        "num_filtered_fragments": len(filtered_fragments),
        "num_representatives": len(representatives),
    }


class Clustering:
    """
    Orchestrator class for representative fragment clustering.

    This class probides a thin wrapper around the clustering workflow,
    primarily to maintain consistency with other COFFEE-PRESC components.
    """

    def __init__(self, settings=None):
        self.__settings: RepclusSettings = settings
        self._stats = {}

    @property
    def settings(self) -> RepclusSettings:
        return self.__settings

    def run(
        self,
    ) -> None:
        """
        Run the representative fragment clustering process.

        Temporary files are created internally for intermediate molecule
        and fragment representations and are cleaned up automatically.

        Parameters
        ----------
        molecules_path : StrPath
            Path to the input molecules SDF file.
        output_representative_path : StrPath
            Path to the output SDF file for representative fragments.
        n_clusters : int
            Number of clusters used for fragment clustering.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            tmp_molecules_path = f"{temp_dir}/tmp_molecules.sdf"
            tmp_fragments_path = f"{temp_dir}/tmp_fragments.sdf"

            molecule_filter: Filter = _build_molecule_filter()
            fragment_filter: Filter = _build_fragment_filter()

            self._stats = _clustering_flow(
                molecules_path=self.__settings.clustering_molecules,
                output_representative_path=self.__settings.clustering_output,
                molecule_filter=molecule_filter,
                fragment_filter=fragment_filter,
                tmp_molecules_path=tmp_molecules_path,
                tmp_fragments_path=tmp_fragments_path,
                n_clusters=self.__settings.n_clusters,
            )
