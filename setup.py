import os
import subprocess

from setuptools import find_packages, setup
from setuptools.command.build_py import build_py


# Custom command to build C++ binaries
class build_cpp_bin(build_py):
    def run(self):
        # Read paths from environment variables
        BOOST_ROOT = os.environ.get("BOOST_ROOT", "")
        OPEN_BABEL_ROOT = os.environ.get("OPEN_BABEL_ROOT", "")
        Debug = os.environ.get("Debug", "N")
        Symbols = os.environ.get("Symbols", "N")
        OpenMP = os.environ.get("OpenMP", "N")
        STATIC = os.environ.get("STATIC", "Y")  # Static linking by default

        CXXFLAGS = ["-std=c++11"]
        if Debug == "Y":
            CXXFLAGS += ["-g", "-O0"]
        else:
            CXXFLAGS += ["-O2"]
        if Symbols == "Y":
            CXXFLAGS += ["-g"]
        if OpenMP == "Y":
            CXXFLAGS += ["-fopenmp"]
        STATICFLAG = ["-s", "-lz", "-static"] if STATIC == "Y" else []

        BOOSTIP = [f"-I{BOOST_ROOT}/include"] if BOOST_ROOT else []
        BOOSTLP = [f"-L{BOOST_ROOT}/lib"] if BOOST_ROOT else []
        OBABELIP = (
            [f"-I{OPEN_BABEL_ROOT}/include/openbabel-2.0"] if OPEN_BABEL_ROOT else []
        )
        OBABELLP = [f"-L{OPEN_BABEL_ROOT}/lib"] if OPEN_BABEL_ROOT else []

        fragquery_src_dir = os.path.join("coffeepresc", "fragquery", "src")
        fragquery_bin_dir = os.path.join("coffeepresc", "fragquery", "bin")
        fragquery_objs_dir = os.path.join(fragquery_src_dir, "objs")

        decompose_src_dir = os.path.join("coffeepresc", "fragment", "decompose", "src")
        decompose_bin_dir = os.path.join("coffeepresc", "fragment", "decompose", "bin")
        decompose_objs_dir = os.path.join(decompose_src_dir, "objs")

        for d in (fragquery_bin_dir, fragquery_objs_dir, decompose_bin_dir, decompose_objs_dir):
            os.makedirs(d, exist_ok=True)

        # Object lists for each executable
        grid_objs = [
            "grid_main.o",
            "main_utils.o",
            "utils.o",
            "infile_reader.o",
            "Molecule.o",
            "Fragment.o",
            "Vector3d.o",
            "Atom.o",
            "AtomInterEnergyGrid.o",
            "InterEnergyGrid.o",
            "EnergyCalculator.o",
            "log_writer_stream.o",
            "OBMol.o",
        ]
        fragquery_objs = [
            "fragquery_main.o",
            "main_utils.o",
            "Vector3d.o",
            "InterEnergyGrid.o",
            "Molecule.o",
            "Fragment.o",
            "EnergyCalculator.o",
            "infile_reader.o",
            "utils.o",
            "AtomInterEnergyGrid.o",
            "FragmentInterEnergyGrid.o",
            "Atom.o",
            "log_writer_stream.o",
            "OBMol.o",
            "QueryGenerator.o",
        ]

        decompose_objs = [
            "decompose_main.o",
            "UnionFindTree.o",
            "utils.o",
            "OBMolReader.o",
            "Converter.o",
            "MoleculeToFragments.o",
            "AtomType.o",
            "Atom.o",
            "Molecule.o",
            "Vector3d.o",
            "Fragment.o",
            "infile_reader.o",
            "log_writer_stream.o",
        ]

        # Build .o files
        def build_obj(obj, src_dir, objs_dir):
            src_name = obj.replace(".o", ".cc")
            src_path = os.path.join(src_dir, src_name)
            obj_path = os.path.join(objs_dir, obj)
            if not os.path.exists(src_path):
                return
            cmd = (
                ["g++"]
                + CXXFLAGS
                + BOOSTIP
                + OBABELIP
                + ["-o", obj_path, "-c", src_path]
            )
            print(f"Compiling {obj}: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)

        for obj in sorted(set(grid_objs + fragquery_objs)):
            build_obj(obj, fragquery_src_dir, fragquery_objs_dir)

        for obj in decompose_objs:
            build_obj(obj, decompose_src_dir, decompose_objs_dir)

        # Link executables
        def link_executable(exe_name, obj_list, objs_dir, bin_dir):
            obj_paths = [os.path.join(objs_dir, o) for o in obj_list]
            exe_path = os.path.join(bin_dir, exe_name)
            cmd = (
                ["g++"]
                + CXXFLAGS
                + ["-o", exe_path]
                + obj_paths
                + BOOSTLP
                + OBABELLP
                + ["-lboost_regex", "-lboost_program_options", "-lopenbabel"]
                + STATICFLAG
            )
            print(f"Linking {exe_name}: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)

        link_executable("atomgrid-gen", grid_objs, fragquery_objs_dir, fragquery_bin_dir)
        link_executable("fragment-query", fragquery_objs, fragquery_objs_dir, fragquery_bin_dir)
        link_executable("decompose", decompose_objs, decompose_objs_dir, decompose_bin_dir)

        super().run()


setup(
    name="coffeepresc",
    version="1.0.0",
    description="COFFEE-PRESC: COmpound Filtering by Fragment pair-based Efficient Evaluation for PRE-SCreening",
    packages=find_packages(),
    install_requires=[
        "numpy>=2.0.0",
        "rdkit==2024.3.5",
        "h5py",
        "pandas",
    ],
    include_package_data=True,
    package_data={
        "coffeepresc.fragquery.bin": ["atomgrid-gen", "fragment-query"],
        "coffeepresc.fragment.decompose.bin": ["decompose"],
    },
    python_requires=">=3.12",
    cmdclass={
        "build_cpp": build_cpp_bin,
    },
    entry_points={
        "console_scripts": [
            "coffeepresc=coffeepresc.main:main",
            "fragquery=coffeepresc.fragquery.cli:main",
            "fbdb=coffeepresc.retrieval.cli:main",
            "cmpdeval=coffeepresc.cmpdeval.cli:main",
            "repclus=coffeepresc.fragment.repclus.cli:main",
            "decompose=coffeepresc.fragment.decompose.cli:main",
        ],
    },
)
