#!/usr/bin/env python3
"""
CHGNet + ASE relax/short-MD CLI with simple provenance and sanity checks.
"""
from __future__ import annotations

import argparse, json, os, time, sys, math, pathlib
from dataclasses import asdict, dataclass

from ase import io
from ase.optimize import BFGS
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units

# Lazy import to avoid heavy cost when --help
def _load_calculator():
    from chgnet.model import CHGNet
    from chgnet.model.dynamics import CHGNetCalculator
    model = CHGNet.load()
    return CHGNetCalculator(model)

@dataclass
class SanitySummary:
    max_force_eVA: float
    energy_eV: float
    md_steps: int
    md_energy_drift_meV_per_atom: float | None
    flagged: bool
    reasons: list[str]

def write_json(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)

def get_versions():
    import chgnet, ase, numpy
    return {"chgnet": chgnet.__version__, "ase": ase.__version__, "numpy": numpy.__version__}

def parse_args():
    p = argparse.ArgumentParser(description="Pre-relax with CHGNet + optional short MD")
    p.add_argument("structure", help="Path to CIF or POSCAR")
    p.add_argument("--outdir", default="runs/out", help="Output directory")
    p.add_argument("--fmax", type=float, default=0.05, help="Relaxation threshold in eV/Å (default 0.05)")
    p.add_argument("--md-steps", type=int, default=0, help="Number of MD steps (0 to skip)")
    p.add_argument("--md-timestep-fs", type=float, default=1.0, help="MD timestep (fs)")
    p.add_argument("--md-temperature-K", type=float, default=300.0, help="MD temperature (K)")
    p.add_argument("--seed", type=int, default=42, help="Random seed for velocities")
    p.add_argument("--drift-threshold-meV-per-atom", type=float, default=2.5, help="Flag if abs(ΔE)/N > threshold over MD")
    return p.parse_args()

def main():
    args = parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    atoms = io.read(args.structure)
    calc = _load_calculator()
    atoms.calc = calc

    # Relax
    opt = BFGS(atoms, logfile=str(outdir / "opt.log"))
    opt.run(fmax=args.fmax)
    energy = float(atoms.get_potential_energy())

    # Save relaxed structure
    io.write(outdir / "relaxed.cif", atoms)
    try:
        io.write(outdir / "relaxed.POSCAR", atoms, format="vasp")
    except Exception:
        pass

    # Sanity: max force
    max_f = float((abs(atoms.get_forces())).max())

    md_drift = None
    if args.md_steps > 0:
        # Initialize velocities, run short MD
        MaxwellBoltzmannDistribution(atoms, temperature_K=args.md_temperature_K, rng=args.seed)
        dyn = Langevin(atoms, args.md_timestep_fs * units.fs, temperature_K=args.md_temperature_K, friction=0.02)
        energies = []
        for step in range(args.md_steps):
            dyn.run(1)
            energies.append(float(atoms.get_potential_energy()))
        # Save trajectory
        io.write(outdir / "md.traj", atoms)
        # Drift per atom
        if len(energies) >= 2:
            md_drift = 1000.0 * (energies[-1] - energies[0]) / len(atoms)  # meV/atom

        # Save simple summary
        write_json({"energies_eV": energies}, outdir / "md_summary.json")

    # Sanity & flags
    flagged = False
    reasons = []
    if max_f > max(0.1, args.fmax * 2):
        flagged = True
        reasons.append(f"High max force after relax: {max_f:.3f} eV/Å")
    if md_drift is not None and abs(md_drift) > args.drift_threshold_meV_per_atom:
        flagged = True
        reasons.append(f"MD energy drift {md_drift:.2f} meV/atom exceeds threshold {args.drift_threshold_meV_per_atom:.2f}")

    sanity = SanitySummary(max_force_eVA=max_f, energy_eV=energy, md_steps=args.md_steps,
                           md_energy_drift_meV_per_atom=md_drift, flagged=flagged, reasons=reasons)
    write_json(asdict(sanity), outdir / "sanity.json")

    # Provenance
    prov = {
        "cli_args": vars(args),
        "versions": get_versions(),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "structure_input": str(args.structure),
        "outputs": {
            "relaxed_cif": str(outdir / "relaxed.cif"),
            "relaxed_poscar": str(outdir / "relaxed.POSCAR"),
            "sanity": str(outdir / "sanity.json"),
        }
    }
    write_json(prov, outdir / "provenance.json")
    print("Done. Outputs in", outdir)

if __name__ == "__main__":
    main()
