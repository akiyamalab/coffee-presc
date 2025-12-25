from __future__ import annotations
from collections import defaultdict
from typing import Self

import numpy as np
from coffeepresc.fragment.repclus.similarity import (
    create_distance_matrix_by_fragment_similarity,
)
from rdkit import Chem
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist


class Clusters:
    """
    Container class representing a set of fragment clusters.

    This class manages multiple `Cluster` instances produced by hierarchical
    clustering and provides utilities for selecting representative fragments
    from each cluster.

    Notes
    -----
    - Clustering is performed externally via the `from_fragments_by_ward` factory.
    - Representative fragments are computed lazily and cached internally.
    """

    clusters: list[Cluster]
    representatives: list[Chem.Mol] | None

    def __init__(self, clusters: list[Cluster]):
        """
        Initialize the Clusters container.

        Parameters
        ----------
        clusters : list[Cluster]
            List of Cluster instances, each representing a fragment group.
        """

        self.clusters = clusters
        self.representatives = None

    @classmethod
    def from_fragments_by_ward(cls, molecules: list[Chem.Mol], n_cluster: int) -> Self:
        """
        Construct clusters from fragment molecules using Ward hierarchical clustering.

        A distance matrix based on fragment similarity is computed, followed by
        Ward linkage and flat clustering with a fixed number of clusters.

        Parameters
        ----------
        molecules : list[Chem.Mol]
            Fragment molecules to be clustered.
        n_cluster : int
            Desired number of clusters.

        Returns
        -------
        Clusters
            A Clusters instance containing clustered fragments.
        """

        # create distance matrix
        distance_matrix: np.ndarray = create_distance_matrix_by_fragment_similarity(
            molecules=molecules
        )
        # ward clustering
        labels: np.ndarray = cls.compute_labels_by_ward(
            distance_matrix=distance_matrix, n_cluster=n_cluster
        )
        # create clusters
        cluster_dict: dict[int, list[Chem.Mol]] = defaultdict(list)
        for mol, label in zip(molecules, labels):
            cluster_dict[label].append(mol)

        # create Cluster instances
        clusters: list[Cluster] = [
            Cluster(molecules=mol_list) for mol_list in cluster_dict.values()
        ]

        return cls(clusters=clusters)

    @staticmethod
    def compute_labels_by_ward(
        distance_matrix: np.ndarray, n_cluster: int
    ) -> np.ndarray:
        """
        Compute cluster labels using Ward hierarchical clustering.

        Parameters
        ----------
        distance_matrix : np.ndarray
            Pairwise distance matrix between fragments.
        n_cluster : int
            Number of clusters to form.

        Returns
        -------
        np.ndarray
            Cluster labels assigned to each fragment.
        """

        distance_array = pdist(distance_matrix)
        linkage_matrix = linkage(distance_array, method="ward")
        labels: np.ndarray = fcluster(linkage_matrix, n_cluster, criterion="maxclust")

        return labels

    def build_representatives(self) -> None:
        """
        Select representative fragments for all clusters.

        The representative of each cluster is defined as the fragment with the
        minimum mean distance to all other fragments in the same cluster.
        """

        representatives: list[Chem.Mol] = [
            cluster.select_representative_fragment() for cluster in self.clusters
        ]

        self.representatives = representatives

    def get_representatives(self) -> list[Chem.Mol]:
        """
        Retrieve representative fragments for each cluster.

        Returns
        -------
        list[Chem.Mol]
            Representative fragments, one per cluster.

        Raises
        ------
        ValueError
            If representatives have not been built yet.
        """

        if self.representatives is None:
            raise ValueError("Representatives have not been built yet.")
        return self.representatives


class Cluster:
    """
    Class representing a single fragment cluster.

    This class stores fragment molecules belonging to the same cluster and
    provides a method to select a representative fragment.
    """

    molecules: list[Chem.Mol]
    representative: Chem.Mol | None

    def __init__(self, molecules: list[Chem.Mol]):
        """
        Class representing a single fragment cluster.

        This class stores fragment molecules belonging to the same cluster and
        provides a method to select a representative fragment.
        """

        self.molecules = molecules
        self.representative = None

    def select_representative_fragment(self) -> Chem.Mol:
        """
        Select a representative fragment for this cluster.

        The representative is chosen as the fragment with the smallest
        average distance to all other fragments in the cluster.

        Returns
        -------
        Chem.Mol
            The selected representative fragment.
        """

        distance_matrix: np.ndarray = create_distance_matrix_by_fragment_similarity(
            molecules=self.molecules
        )

        mean_distances: np.ndarray = np.mean(distance_matrix, axis=1)
        representative_index: int = np.argmin(mean_distances)
        self.representative = self.molecules[representative_index]

        return self.representative


if __name__ == "__main__":
    mols = list(Chem.SDMolSupplier("example/fragments.sdf"))
    clusters = Clusters.from_fragments_by_ward(mols, n_cluster=5)
    clusters.build_representatives()
    reps = clusters.get_representatives()
    for i, rep in enumerate(reps):
        print(f"Representative of cluster {i+1}: {Chem.MolToSmiles(rep)}")
