# COFFEE-PRESC
COFFEE-PRESC (COmpound Filtering by Fragment pair-based Efficient Evaluation for PRE-SCreening): A fast pre-screening method using compound retrieval by pairwise positional relationship of representative fragments

## Highlights
- Fragment docking and query enumeration backed by native C++ tools
- Fragment clustering preprocessing via `repclus` to generate representative fragments
- Fragment decomposition via `decompose`
- Efficient similarity-based retrieval on an HDF5 store
- End-to-end CLI (`coffeepresc`) plus modular CLIs (`fragquery`, `fbdb`, `cmpdeval`, `repclus`, `decompose`)
- Configurable via a single TOML settings file

## Requirements
- OS: Linux (tested on Ubuntu 20.04)
- Python: 3.12+
- Python packages: `numpy>=2.0.0`, `rdkit==2024.3.5`, `h5py`, `pandas`, `scipy`
- Optional/for rebuilding C++ binaries:
  - `g++` with C++11
  - Boost (regex, program_options)
  - Open Babel 2.4.1
    - **Note: Open Babel 3.x is not supported**

### Recommended: Use Dev Container / Docker
This repository includes ready-to-use dev containers (`.devcontainer/` with Dockerfile + devcontainer.json). **We recommend using Dev Containers or Docker** to avoid manual dependency setup and ensure consistency across environments.

## Installation
```bash
# From the repository root
python -m pip install .
# or for development
python -m pip install -e .
```

Required for `fragquery`, `coffeepresc`, `repclus`, and `decompose`: build native C++ tools
```bash
# Builds atomgrid-gen and fragment-query into coffeepresc/fragquery/bin
python setup.py build_cpp
```
Without this step, `fragquery` and the full pipeline `coffeepresc` will not run.
Note: Building requires Boost and Open Babel 2.4.1. If installed in non-standard locations, set `BOOST_ROOT` and `OPEN_BABEL_ROOT`.

## Quick Start
```bash
# 1) Preprocessing (database build)
fbdb create -s example/example.toml --conformers example/conformers.sdf.gz --log example/fbdb_create.log

# 2) Run the end-to-end pipeline (fragquery → fbdb search → cmpdeval)
coffeepresc -s example/example.toml

# Or run step-by-step
fragquery -s example/example.toml
fbdb search -s example/example.toml
cmpdeval -s example/example.toml
```

## Configuration (TOML)
All CLIs accept `-s/--setting` pointing to a TOML file. Configuration precedence: **command-line arguments > TOML > defaults**. Typical workflow: set base configuration in TOML, override or add specific parameters via command-line flags as needed.

Example:
```toml
# COFFEE-PRESC files
## Input
receptor = "example/receptor.pdb"
fragments = "example/fragments.sdf"
## Output
output = "example/output.csv"
## Intermediate (optional in `coffeepresc`)
query = "example/query.csv"
matched = "example/matched.csv"

# Logger option
log = "example/coffeepresc.log"

# For fragquery
docking_config = "example/docking_box.conf"

# For fbdb
storage = "example/db.h5"
similarity_th = 0.45

# For cmpdeval
penalty_coef = 8.0
```
Notes:
- If `query`/`matched`/`grid` are omitted, temporary files/folders will be created.
- `storage` must exist before calling `fbdb search` or the `coffeepresc` pipeline.

### Docking Configuration File Format
The `docking_config` file (e.g., `example/docking_box.conf`) defines the docking region using a simple key-value format:
```
INNERBOX 16, 16, 16
OUTERBOX 26, 26, 26
BOX_CENTER 46.9322, -19.9701, 102.5602
SCORING_PITCH 0.25, 0.25, 0.25
```
- `INNERBOX`: Fragment grid dimensions (Å) for fragment placement (x, y, z)
- `OUTERBOX`: Atom grid dimensions (Å) defining the extended search region (x, y, z)
- `BOX_CENTER`: Center coordinates (Å) of the docking box (x, y, z)
- `SCORING_PITCH`: Grid spacing (Å) for energy calculation (x, y, z)

## Command-Line Interfaces
Six console scripts are provided via `setup.py`:

### coffeepresc
Runs the full pipeline: fragment grid + query enumeration → database search → scoring.
Requires native C++ tools built via `python setup.py build_cpp`.
```bash
coffeepresc -s example/example.toml
```
For available options (TOML parameters and command-line overrides), see `coffeepresc --help`.

### fragquery
Generates grids and enumerates fragment queries using native C++ tools.
Requires native C++ tools built via `python setup.py build_cpp`.
```bash
fragquery -s example/example.toml
```
For available options (TOML parameters and command-line overrides), see `fragquery --help`.

### fbdb
Manages the fragment-based retrieval database.
```bash
# Create database from conformers and representative fragments
fbdb create -s example/example.toml --conformers example/conformers.sdf.gz

# Search matched compounds for a query
fbdb search -s example/example.toml
```
For available options (TOML parameters and command-line overrides), see `fbdb create --help` or `fbdb search --help`.

### cmpdeval
Scores matched compounds and writes final results.
```bash
cmpdeval -s example/example.toml
```
For available options (TOML parameters and command-line overrides), see `cmpdeval --help`.

### repclus
Representative fragment clustering (preprocessing; produces a representative fragments SDF for downstream steps such as `fragquery` and `fbdb create` via `--fragments`). Requires native C++ tools built via `python setup.py build_cpp`. Does not support `-s/--setting` TOML; pass arguments explicitly. For details, see `repclus --help`.

### decompose
Fragment decomposition (preprocessing; produces fragments SDF from molecules for `fbdb create`). Requires native C++ tools built via `python setup.py build_cpp`. Does not support `-s/--setting` TOML; pass arguments explicitly. For details, see `decompose --help`.
- **Note**: The output from `decompose` requires conformer generation using an external tool (e.g., Omega). The resulting conformers file will be used later as input to `fbdb create` (via `--conformers`).

## Reference
Shimizu M, Yoneyama S, Yanagisawa K, Akiyama Y. [COFFEE-PRESC: A Fast Prescreening Method Using Compound Retrieval by Pairwise Positional Relationship of Representative Fragments](https://pubs.acs.org/doi/10.1021/acs.jcim.5c03067), *Journal of Chemical Information and Modeling*, 66(8):4672–4684, 2026. doi: 10.1021/acs.jcim.5c03067
