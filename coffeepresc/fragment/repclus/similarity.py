import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdMolDescriptors

MQN_NUM = 42
MORGAN_RADIUS = 2
MORGAN_NBITS = 1024


def create_distance_matrix_by_molecular_quantum_numbers_similarity(
    molecules: list[Chem.Mol],
) -> np.ndarray:
    """
    Compute a distance matrix based on Molecular Quantum Numbers (MQNs).

    For each molecule, a 42-dimensional MQN descriptor is computed.
    Pairwise similarities are calculated as an inverse function of the
    L1 distance between MQN vectors, normalized by the number of MQN
    features. The final output is converted to a distance matrix.

    For small molecule sets, a fully vectorized implementation is used
    for efficiency. For large sets, a nested-loop fallback is applied
    to reduce memory consumption.

    Parameters
    ----------
    molecules : list[Chem.Mol]
        List of RDKit molecule objects.

    Returns
    -------
    np.ndarray
        A square distance matrix of shape (N, N), where N is the number
        of molecules. Each entry ranges from 0 (identical) to 1
        (maximally dissimilar).
    """

    mqns = np.array([rdMolDescriptors.MQNs_(molecule) for molecule in molecules])

    if mqns.shape[0] < 1000:  # for memory efficiency
        similarity_matrix = np.abs(mqns[:, None] - mqns)
        similarity_matrix = 1 / (1 + np.sum(similarity_matrix, axis=2) / MQN_NUM)
    else:
        similarity_matrix = np.array(
            [
                [1 / (1 + np.sum(np.abs(mqns_i - mqns_j)) / MQN_NUM) for mqns_j in mqns]
                for mqns_i in mqns
            ]
        )

    return 1 - similarity_matrix


def create_distance_matrix_by_morgan_fingerprint_similarity(
    molecules: list[Chem.Mol],
) -> np.ndarray:
    """
    Compute a distance matrix based on Morgan fingerprint similarity.

    Each molecule is encoded as a Morgan (ECFP-like) fingerprint with
    a fixed radius and bit length. Pairwise similarities are computed
    using the Tanimoto coefficient and converted into distances.

    Distance is defined as:
        distance = 1 - TanimotoSimilarity

    Parameters
    ----------
    molecules : list[Chem.Mol]
        List of RDKit molecule objects.

    Returns
    -------
    np.ndarray
        Square distance matrix of shape (N, N), where N is the number
        of molecules. Values range from 0 (identical fingerprints) to
        1 (no shared features).
    """

    morgan_fp = [
        AllChem.GetMorganFingerprintAsBitVect(
            molecule, MORGAN_RADIUS, nBits=MORGAN_NBITS
        )
        for molecule in molecules
    ]
    distance_matrix = 1 - np.array(
        [
            DataStructs.BulkTanimotoSimilarity(morgan_fp[i], morgan_fp)
            for i in range(len(morgan_fp))
        ]
    )

    return distance_matrix


def create_distance_matrix_by_fragment_similarity(
    molecules: list[Chem.Mol],
) -> np.ndarray:
    """
    Compute a fragment-level distance matrix by combining MQN and
    Morgan fingerprint similarities.

    Two complementary distance matrices are computed:
    - MQN-based distance (global compositional similarity)
    - Morgan fingerprint distance (local substructure similarity)

    The final distance matrix is obtained by averaging the two:
        distance = (mqn_distance + morgan_distance) / 2

    This hybrid metric balances global physicochemical features and
    local structural patterns, making it suitable for fragment
    clustering and representative selection.

    Parameters
    ----------
    molecules : list[Chem.Mol]
        List of RDKit molecule objects.

    Returns
    -------
    np.ndarray
        Square distance matrix of shape (N, N), where N is the number
        of fragments. Values range from 0 (highly similar) to 1
        (highly dissimilar).
    """

    mqn_distance_matrix = (
        create_distance_matrix_by_molecular_quantum_numbers_similarity(molecules)
    )
    morgan_distance_matrix = create_distance_matrix_by_morgan_fingerprint_similarity(
        molecules
    )

    return (mqn_distance_matrix + morgan_distance_matrix) / 2


if __name__ == "__main__":
    mols = list(Chem.SDMolSupplier("example/fragments.sdf"))
    distance_matrix = create_distance_matrix_by_fragment_similarity(mols)
    print(distance_matrix)
