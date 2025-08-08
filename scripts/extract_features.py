#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array


def parse_lipid_list(csv_path: str):
    df = pd.read_csv(csv_path)
    lipids = sorted(df['lipid'].unique().tolist())
    return lipids


def compute_midplane(u: mda.Universe, lipids):
    # Wähle Head- oder alle Atome der Membran-Lipide und berechne Midplane als Median der z-Koordinate
    sel = " or ".join([f"resname {lip}" for lip in lipids])
    mem = u.select_atoms(sel)
    if mem.n_atoms == 0:
        raise ValueError("Keine Membran-Atome gefunden; prüfe Lipid-Resnamen und Topologie.")
    zvals = mem.positions[:, 2]
    z_mid = np.median(zvals)
    return z_mid


def count_contacts(u: mda.Universe, protein_sel: str, membrane_sel: str, cutoff=6.0):
    # cutoff in Angström (MDAnalysis default)
    prot = u.select_atoms(protein_sel)
    mem = u.select_atoms(membrane_sel)
    if prot.n_atoms == 0 or mem.n_atoms == 0:
        return 0
    d = distance_array(prot.positions, mem.positions)
    return int((d < cutoff).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tpr', required=True)
    ap.add_argument('--xtc', required=True)
    ap.add_argument('--lipid_csv', required=True)
    ap.add_argument('--out_csv', required=True)
    ap.add_argument('--last_ns', type=float, default=5.0)
    args = ap.parse_args()

    u = mda.Universe(args.tpr, args.xtc)
    dt_ps = float(u.trajectory.dt) if hasattr(u.trajectory, 'dt') else 2.0
    total_ps = dt_ps * len(u.trajectory)
    start_ps = max(0.0, total_ps - args.last_ns * 1000.0)

    lipids = parse_lipid_list(args.lipid_csv)

    # Selektoren
    protein_sel = 'protein'
    membrane_sel = " or ".join([f"resname {lip}" for lip in lipids])

    rows = []
    for ts in u.trajectory:
        if ts.time < start_ps:
            continue
        z_mid = compute_midplane(u, lipids)
        prot = u.select_atoms(protein_sel)
        z_prot = prot.center_of_mass()[2]
        insertion_nm = abs((z_prot - z_mid) / 10.0)  # Angström -> nm
        contacts = count_contacts(u, protein_sel, membrane_sel, cutoff=6.0)
        rows.append({
            'time_ps': ts.time,
            'insertion_nm': insertion_nm,
            'contacts_6A': contacts
        })

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    print(f"[OK] Features geschrieben: {args.out_csv}")


if __name__ == '__main__':
    main()
