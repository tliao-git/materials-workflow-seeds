# chgnet_relax_cli

A minimal command‑line tool to pre‑relax a crystal structure with CHGNet + ASE, optionally run a short MD "shake test," and write simple provenance & sanity checks.

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install --upgrade pip
pip install chgnet>=0.3.4 ase>=3.22.1 numpy>=1.23 matplotlib>=3.7
# CHGNet depends on PyTorch; if needed:
# pip install torch --index-url https://download.pytorch.org/whl/cpu
```

## Usage

```bash
python relax.py structures/Fe.cif --md-steps 100 --fmax 0.05 --outdir runs/fe_test
```

Outputs:
- `relaxed.cif` and `relaxed.POSCAR`
- `md.traj` (if MD requested) and `md_summary.json`
- `provenance.json` with package versions, CLI args, and timestamps
- `sanity.json` with simple checks (max force, energy drift)

## Notes
- Input can be CIF or POSCAR. Multi‑frame inputs will use the first frame.
- Sanity thresholds are conservative and can be tuned with CLI flags.
