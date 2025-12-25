from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


class Filter:
    """
    Configurable molecule filter for preprocessing RDKit molecules.

    This class provides a lightweight mechanism to compose and apply
    multiple molecule-level filtering rules. Each filter is registered
    as a callable and applied sequentially to a list of RDKit `Chem.Mol`
    objects.

    The filter is primarily used to remove chemically undesirable or
    out-of-scope molecules and fragments prior to representative
    fragment generation.

    Notes
    -----
    - Filters are applied in the order they are added.
    - The class is intentionally simple and stateful, prioritizing
      readability over extensibility.
    """

    ALLOWED_ELEMENTS = [
        "H",
        "C",
        "N",
        "O",
        "F",
        "P",
        "S",
        "Cl",
        "Br",
        "I",
    ]
    verbose: bool = False

    def __init__(self):
        self.filters: list[tuple[str, callable[[list[Chem.Mol]], list[Chem.Mol]]]] = []

    def add_remove_Hs(self) -> None:
        """
        Remove hydrogen atoms from molecules.

        Parameters
        ----------
        allowed_elements : list[str], optional
            List of allowed atomic symbols (default: common organic elements).

        Notes
        -----
        - Molecules containing any disallowed element are removed.
        - `None` molecules are ignored.
        """

        def remove_Hs_from_molecules(molecules: list[Chem.Mol]) -> list[Chem.Mol]:
            molecules = [molecule for molecule in molecules if molecule is not None]
            filtered_molecules = [Chem.RemoveHs(molecule) for molecule in molecules]
            if self.verbose:
                print(
                    f"Removed Hs from {len(molecules)} molecules."
                )
            return filtered_molecules

        self.filters.append(("removing Hs", remove_Hs_from_molecules))

    def add_filter_by_allowed_elements(
        self, allowed_elements: list[str] = ALLOWED_ELEMENTS
    ) -> None:
        """
        Filter molecules by allowed elements.

        Parameters
        ----------
        allowed_elements : list[str], optional
            List of allowed atomic symbols (default: common organic elements).

        Notes
        -----
        - Molecules containing any disallowed element are removed.
        - `None` molecules are ignored.
        """

        def filter_molecules_by_allowed_elements(
            molecules: list[Chem.Mol],
        ) -> list[Chem.Mol]:
            molecules = [molecule for molecule in molecules if molecule is not None]
            filtered_molecules = [
                molecule
                for molecule in molecules
                if all(
                    atom.GetSymbol() in allowed_elements for atom in molecule.GetAtoms()
                )
            ]
            if self.verbose:
                print(
                    f"filtering allowed elements: {len(molecules)} -> {len(filtered_molecules)}"
                )
            return filtered_molecules

        self.filters.append(
            ("filtering by allowed elements", filter_molecules_by_allowed_elements)
        )

    def add_filter_by_isotope_atoms(self) -> None:
        """
        Filter molecules containing isotope atoms.

        Parameters
        ----------
        allowed_elements : list[str], optional
            List of allowed atomic symbols (default: common organic elements).

        Notes
        -----
        - Molecules containing any disallowed element are removed.
        - `None` molecules are ignored.
        """

        def filter_molecules_by_isotope_atoms(
            molecules: list[Chem.Mol],
        ) -> list[Chem.Mol]:
            molecules = [molecule for molecule in molecules if molecule is not None]
            filtered_molecules = [
                molecule
                for molecule in molecules
                if all(atom.GetIsotope() == 0 for atom in molecule.GetAtoms())
            ]
            if self.verbose:
                print(
                    f"filtering by isotope atoms: {len(molecules)} -> {len(filtered_molecules)}"
                )
            return filtered_molecules

        self.filters.append(
            ("filtering by isotope atoms", filter_molecules_by_isotope_atoms)
        )

    def add_filter_by_ring_size(self, ring_size: int) -> None:
        """
        Filter molecules by ring size.

        Parameters
        ----------
        allowed_elements : list[str], optional
            List of allowed atomic symbols (default: common organic elements).

        Notes
        -----
        - Molecules containing any disallowed element are removed.
        - `None` molecules are ignored.
        """

        def filter_molecules_by_ring_size(molecules: list[Chem.Mol]) -> list[Chem.Mol]:
            molecules = [molecule for molecule in molecules if molecule is not None]
            filtered_molecules = [
                molecule
                for molecule in molecules
                if all(
                    molecule.GetRingInfo().MinAtomRingSize(atom.GetIdx()) < ring_size
                    for atom in molecule.GetAtoms()
                )
            ]
            if self.verbose:
                print(
                    f"filtering by ring size: {len(molecules)} -> {len(filtered_molecules)}"
                )
            return filtered_molecules

        self.filters.append(("filtering by ring size", filter_molecules_by_ring_size))

    def add_filter_by_ring_number(self, ring_number: int) -> None:
        """
        Filter molecules by ring number.

        Parameters
        ----------
        allowed_elements : list[str], optional
            List of allowed atomic symbols (default: common organic elements).

        Notes
        -----
        - Molecules containing any disallowed element are removed.
        - `None` molecules are ignored.
        """

        def filter_molecules_by_ring_number(
            molecules: list[Chem.Mol],
        ) -> list[Chem.Mol]:
            molecules = [molecule for molecule in molecules if molecule is not None]
            filtered_molecules = [
                molecule
                for molecule in molecules
                if molecule.GetRingInfo().NumRings() < ring_number
            ]
            if self.verbose:
                print(
                    f"filtering by ring number: {len(molecules)} -> {len(filtered_molecules)}"
                )
            return filtered_molecules

        self.filters.append(
            ("filtering by ring number", filter_molecules_by_ring_number)
        )

    def add_filter_by_aromatic_ring_number(self, ring_number: int) -> None:
        """
        Filter molecules by aromatic ring number.

        Parameters
        ----------
        allowed_elements : list[str], optional
            List of allowed atomic symbols (default: common organic elements).

        Notes
        -----
        - Molecules containing any disallowed element are removed.
        - `None` molecules are ignored.
        """

        def filter_molecules_by_aromatic_ring_number(
            molecules: list[Chem.Mol],
        ) -> list[Chem.Mol]:
            molecules = [molecule for molecule in molecules if molecule is not None]
            filtered_molecules = [
                molecule
                for molecule in molecules
                if rdMolDescriptors.CalcNumAromaticRings(molecule) < ring_number
            ]
            if self.verbose:
                print(
                    f"filtering by aromatic ring number: {len(molecules)} -> {len(filtered_molecules)}"
                )
            return filtered_molecules

        self.filters.append(
            (
                "filtering by aromatic ring number",
                filter_molecules_by_aromatic_ring_number,
            )
        )

    def add_filter_by_non_aromatic_ring_number(self, ring_number: int) -> None:
        """
        Filter molecules by non-aromatic ring number.

        Parameters
        ----------
        allowed_elements : list[str], optional
            List of allowed atomic symbols (default: common organic elements).

        Notes
        -----
        - Molecules containing any disallowed element are removed.
        - `None` molecules are ignored.
        """

        def filter_molecules_by_non_aromatic_ring_number(
            molecules: list[Chem.Mol],
        ) -> list[Chem.Mol]:
            molecules = [molecule for molecule in molecules if molecule is not None]
            filtered_molecules = [
                molecule
                for molecule in molecules
                if rdMolDescriptors.CalcNumAliphaticRings(molecule) < ring_number
            ]
            if self.verbose:
                print(
                    f"filtering by non-aromatic ring number: {len(molecules)} -> {len(filtered_molecules)}"
                )
            return filtered_molecules

        self.filters.append(
            (
                "filtering by non-aromatic ring number",
                filter_molecules_by_non_aromatic_ring_number,
            )
        )

    def add_filter_by_heavy_atoms(self, atom_number: int) -> None:
        """
        Filter molecules ( which have no rings ) by heavy atoms.

        Parameters
        ----------
        allowed_elements : list[str], optional
            List of allowed atomic symbols (default: common organic elements).

        Notes
        -----
        - Molecules containing any disallowed element are removed.
        - `None` molecules are ignored.
        """

        def filter_molecules_by_heavy_atoms(
            molecules: list[Chem.Mol],
        ) -> list[Chem.Mol]:
            molecules = [molecule for molecule in molecules if molecule is not None]
            filtered_molecules = [
                molecule
                for molecule in molecules
                if molecule.GetNumHeavyAtoms() < atom_number
                or molecule.GetRingInfo().NumRings() > 0
            ]
            if self.verbose:
                print(
                    f"filtering by heavy atoms: {len(molecules)} -> {len(filtered_molecules)}"
                )
            return filtered_molecules

        self.filters.append(
            ("filtering by heavy atoms", filter_molecules_by_heavy_atoms)
        )

    def add_filter_single_heavy_atom_molecules(self) -> None:
        """
        Filter molecules with single heavy atom.

        Parameters
        ----------
        allowed_elements : list[str], optional
            List of allowed atomic symbols (default: common organic elements).

        Notes
        -----
        - Molecules containing any disallowed element are removed.
        - `None` molecules are ignored.
        """

        def filter_single_heavy_atom_molecules(
            molecules: list[Chem.Mol],
        ) -> list[Chem.Mol]:
            molecules = [molecule for molecule in molecules if molecule is not None]
            filtered_molecules = [
                molecule for molecule in molecules if molecule.GetNumHeavyAtoms() != 1
            ]
            if self.verbose:
                print(
                    f"filtering single heavy atom molecules: {len(molecules)} -> {len(filtered_molecules)}"
                )
            return filtered_molecules

        self.filters.append(
            (
                "filtering single heavy atom molecules",
                filter_single_heavy_atom_molecules,
            )
        )

    def add_filter_molecules_without_carbon(self) -> None:
        """
        Filter molecules without carbon atoms.

        Parameters
        ----------
        allowed_elements : list[str], optional
            List of allowed atomic symbols (default: common organic elements).

        Notes
        -----
        - Molecules containing any disallowed element are removed.
        - `None` molecules are ignored.
        """

        def filter_molecules_without_carbon(
            molecules: list[Chem.Mol],
        ) -> list[Chem.Mol]:
            molecules = [molecule for molecule in molecules if molecule is not None]
            filtered_molecules = [
                molecule
                for molecule in molecules
                if any(atom.GetSymbol() == "C" for atom in molecule.GetAtoms())
            ]
            if self.verbose:
                print(
                    f"filtering by carbon atoms: {len(molecules)} -> {len(filtered_molecules)}"
                )
            return filtered_molecules

        self.filters.append(
            ("filtering molecules without carbon", filter_molecules_without_carbon)
        )

    def apply(self, molecules: list[Chem.Mol]) -> list[Chem.Mol]:
        """
        Apply all registered filters sequentially to a list of molecules.

        Each filter function stored in this instance is applied in the
        order it was added. The output of one filter becomes the input
        to the next filter.

        Parameters
        ----------
        molecules : list[Chem.Mol]
            List of RDKit molecules to be filtered.

        Returns
        -------
        list[Chem.Mol]
            Molecules that remain after all filters have been applied.
        """

        for _, filter_function in self.filters:
            molecules = filter_function(molecules)
        return molecules

    def clear(self) -> None:
        """
        Remove all registered filters.

        After calling this method, the filter pipeline becomes empty
        and `apply()` will return the input molecules unchanged.
        """

        self.filters = []

    def print_filters(self) -> None:
        """
        Print the names of all registered filters in order.

        This method is intended for debugging and inspection, allowing
        users to confirm which filters are currently active and in
        what order they will be applied.
        """

        for filter_name, _ in self.filters:
            print(f"- {filter_name}")

    def set_fragment_filters(self) -> None:
        """
        Register the default set of filters for fragment preprocessing.

        These filters are designed for representative fragment generation
        and enforce constraints such as:
        - allowed elements only,
        - reasonable ring sizes and ring counts,
        - exclusion of isotopes,
        - removal of trivial or chemically uninformative fragments.

        The filters are appended in a predefined order.
        """

        self.add_remove_Hs()
        self.add_filter_by_allowed_elements()
        self.add_filter_by_ring_size(ring_size=7)
        self.add_filter_by_aromatic_ring_number(ring_number=4)
        self.add_filter_by_non_aromatic_ring_number(ring_number=3)
        self.add_filter_by_ring_number(ring_number=4)
        self.add_filter_by_heavy_atoms(atom_number=7)
        self.add_filter_single_heavy_atom_molecules()
        self.add_filter_by_isotope_atoms()

    def set_ligand_filters(self) -> None:
        """
        Register the default set of filters for ligand-level preprocessing.

        These filters remove molecules that are chemically invalid or
        unsuitable as ligands, such as:
        - molecules containing disallowed elements,
        - molecules with isotope atoms,
        - molecules without carbon atoms.

        The filters are appended in a predefined order.
        """

        self.add_filter_by_allowed_elements()
        self.add_filter_by_isotope_atoms()
        self.add_filter_molecules_without_carbon()


if __name__ == "__main__":
    filter = Filter()
    filter.add_remove_Hs()
    filter.add_filter_by_allowed_elements()
    filter.add_filter_by_ring_size(ring_size=7)
    filter.add_filter_by_aromatic_ring_number(ring_number=4)
    filter.add_filter_by_non_aromatic_ring_number(ring_number=3)
    filter.add_filter_by_ring_number(ring_number=4)
    filter.add_filter_by_heavy_atoms(atom_number=7)
    filter.add_filter_single_heavy_atom_molecules()
    filter.add_filter_by_isotope_atoms()
    filter.print_filters()

    mols = list(Chem.SDMolSupplier("example/fragments.sdf"))
    filtered_mols = filter.apply(mols)
    print(f"Input molecules: {len(mols)}")
    print(f"Filtered molecules: {len(filtered_mols)}")
