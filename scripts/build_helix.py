#!/usr/bin/env python3
import sys
import os
import argparse

try:
    import PeptideBuilder
    from Bio.PDB import PDBIO
    from Bio.PDB.Polypeptide import PPBuilder
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Atom import Atom
    from Bio.PDB.PDBParser import PDBParser
except Exception as e:
    sys.stderr.write("[ERROR] Benötigt 'PeptideBuilder' und 'biopython'. Bitte in der Umgebung installieren.\n")
    raise

# Standard-Helix-Winkel für alpha-Helix
PHI = -57.8
PSI_IM1 = -47.0

AMINO3 = {
    'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU','Q':'GLN','G':'GLY','H':'HIS',
    'I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'
}

def build_alpha_helix(sequence: str, out_pdb: str):
    seq3 = [AMINO3[aa] for aa in sequence]
    struct = PeptideBuilder.initialize_res(seq3[0])
    for aa in seq3[1:]:
        PeptideBuilder.add_residue(struct, aa, PHI, PSI_IM1)
    # Schreibe PDB
    io = PDBIO()
    io.set_structure(struct)
    io.save(out_pdb)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sequence", required=True, help="Peptidsequenz, 1-Letter")
    ap.add_argument("--out", required=True, help="Ziel-PDB-Datei")
    args = ap.parse_args()

    seq = args.sequence.strip().upper()
    out_pdb = args.out
    os.makedirs(os.path.dirname(out_pdb), exist_ok=True)
    build_alpha_helix(seq, out_pdb)
    print(f"[OK] Helix-PDB geschrieben: {out_pdb}")
