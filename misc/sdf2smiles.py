import os
import glob
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen

def sdf_folder_to_csv(folder_path, output_csv):
    records = []
    sdf_files = glob.glob(os.path.join(folder_path, "*.sdf"))

    for sdf_file in sdf_files:
        suppl = Chem.SDMolSupplier(sdf_file, sanitize=True, removeHs=False)
        for mol in suppl:
            if mol is None:
                continue
            smiles = Chem.MolToSmiles(mol, canonical=True)
            inchi_key = Chem.MolToInchiKey(mol)
            mol_weight = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            heavy_atoms = mol.GetNumHeavyAtoms()
            records.append({
                "SMILES": smiles,
                "InChIKey": inchi_key,
                "MolWeight": mol_weight,
                "LogP": logp,
                "HeavyAtoms": heavy_atoms,
                "SourceFile": os.path.basename(sdf_file)
            })

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    print(f"Processed {len(records)} molecules from {len(sdf_files)} SDF files. Output written to {output_csv}.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python sdf2smiles.py input_folder output.csv")
        sys.exit(1)
    input_folder = sys.argv[1]
    output_csv = sys.argv[2]
    sdf_folder_to_csv(input_folder, output_csv)

