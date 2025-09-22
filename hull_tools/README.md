# hull_tools

Compute formation energies and convex-hull stability from a simple CSV, using pymatgen. Saves a results CSV and an optional hull plot (binary/ternary only).

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install --upgrade pip
pip install pymatgen>=2024.2.8 numpy>=1.23 pandas>=2.0 matplotlib>=3.7
```

## Input CSV

Required columns:
- `composition` (e.g., Fe2O3, Li1.0Co1.0O2)
- `energy_per_atom_eV` (DFT energy per atom, eV)

Optional:
- `label` (any string to identify a structure)

## Usage

```bash
python hull.py examples/energies.csv --plot hull.png --out results.csv
```

Outputs:
- results CSV with formation energy (eV/atom) and distance to convex hull (eV/atom)
- an optional hull plot (for binaries/ternaries)
