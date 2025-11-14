import sys
import os
import glob
import pandas as pd

try:
    from openbabel import pybel
except ImportError:
    import pybel  # for older Open Babel/Pybel installs

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def pdbqt_folder_to_csv(input_folder, output_csv):
    records = []
    pdbqt_files = glob.glob(os.path.join(input_folder, "*.pdbqt"))
    for pdbqt_file in pdbqt_files:
        for mol in pybel.readfile("pdbqt", pdbqt_file):
            # Convert PDBQT to canonical SMILES using Open Babel
            smiles = mol.write("can").strip()
            # The rest of the properties using RDKit
            try:
                rd_mol = Chem.MolFromSmiles(smiles)
                if rd_mol:
                    inchi_key = Chem.MolToInchiKey(rd_mol)
                    mol_weight = Descriptors.MolWt(rd_mol)
                    logp = Crippen.MolLogP(rd_mol)
                    heavy_atoms = rd_mol.GetNumHeavyAtoms()
                else:
                    inchi_key = ''
                    mol_weight = ''
                    logp = ''
                    heavy_atoms = ''
            except Exception as e:
                inchi_key = ''
                mol_weight = ''
                logp = ''
                heavy_atoms = ''
            records.append({
                "SMILES": smiles,
                "InChIKey": inchi_key,
                "MolWeight": mol_weight,
                "LogP": logp,
                "HeavyAtoms": heavy_atoms,
                "SourceFile": os.path.basename(pdbqt_file)
            })
    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    print(f"Processed {len(records)} molecules from {len(pdbqt_files)} PDBQT files. Output written to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pdbqt2smiles_rdkit.py input_folder output.csv")
        sys.exit(1)
    input_folder = sys.argv[1]
    output_csv = sys.argv[2]
    pdbqt_folder_to_csv(input_folder, output_csv)

