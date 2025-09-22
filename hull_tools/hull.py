#!/usr/bin/env python3
"""
Convex-hull analysis from CSV using pymatgen.
"""
from __future__ import annotations

import argparse, pandas as pd, numpy as np
from pymatgen.core import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.analysis.phase_diagram import PDPlotter
import matplotlib.pyplot as plt

def load_entries(csv_path: str):
    df = pd.read_csv(csv_path)
    if "composition" not in df.columns or "energy_per_atom_eV" not in df.columns:
        raise ValueError("CSV must have columns: composition, energy_per_atom_eV")
    entries = []
    for _, row in df.iterrows():
        comp = Composition(str(row["composition"]))
        epa = float(row["energy_per_atom_eV"])
        entries.append(PDEntry(comp, epa))
    return df, entries

def analyze(csv_path: str, out_csv: str, plot_path: str | None):
    df, entries = load_entries(csv_path)
    pdgm = PhaseDiagram(entries)
    # Compute formation energies and distances
    fenergies = []
    dists = []
    for _, row in df.iterrows():
        comp = Composition(str(row["composition"]))
        entry = PDEntry(comp, float(row["energy_per_atom_eV"]))
        fenergies.append(pdgm.get_form_energy_per_atom(entry))
        dists.append(pdgm.get_e_above_hull(entry))
    df_out = df.copy()
    df_out["formation_energy_per_atom_eV"] = fenergies
    df_out["distance_to_hull_eV"] = dists
    df_out.to_csv(out_csv, index=False)
    print("Wrote", out_csv)

    if plot_path:
        n_elts = len(pdgm.elements)
        if n_elts == 2 or n_elts == 3:
            plotter = PDPlotter(pdgm, show_unstable=True)
            plotter.write_image(plot_path)
            print("Wrote plot", plot_path)
        else:
            print("Plotting supported only for binary/ternary systems. Skipping.")

def main():
    ap = argparse.ArgumentParser(description="Convex-hull analysis with pymatgen")
    ap.add_argument("csv", help="Input CSV with columns: composition, energy_per_atom_eV, [label]")
    ap.add_argument("--out", default="hull_results.csv", help="Output CSV path")
    ap.add_argument("--plot", default=None, help="Optional image path for hull plot (binary/ternary only)")
    args = ap.parse_args()
    analyze(args.csv, args.out, args.plot)

if __name__ == "__main__":
    main()
