import sys
from rdkit import Chem
import pandas as pd

def sdf_to_csv(input_sdf_path, output_csv_path):
    """
    Converts an SDF file to a CSV containing
    MolecularID, SMILES, and HeavyAtoms columns.

    Args:
        input_sdf_path (str): Path to the input SDF file.
        output_csv_path (str): Path to the output CSV file.
    """
    mol_ids = []
    smiles_list = []
    heavy_atoms_list = []

    with open(input_sdf_path, 'r') as sdf_file:
        lines = sdf_file.readlines()

    i = 0
    n_lines = len(lines)
    while i < n_lines:
        # Get molecular ID from the first line
        if lines[i].strip() == '':
            i += 1
            continue
        mol_id = lines[i].strip()
        mol_block = [mol_id + '\n']

        # Read block until terminator "$$$$"
        i += 1
        while i < n_lines and lines[i].strip() != '$$$$':
            mol_block.append(lines[i])
            i += 1

        # Skip terminator lines
        while i < n_lines and lines[i].strip() == '$$$$':
            i += 1

        mol_str = ''.join(mol_block)
        mol = Chem.MolFromMolBlock(mol_str, sanitize=True, removeHs=False)
        if mol:
            smiles = Chem.MolToSmiles(mol)
            heavy_atoms = mol.GetNumHeavyAtoms()
        else:
            smiles = ''
            heavy_atoms = 0

        mol_ids.append(mol_id)
        smiles_list.append(smiles)
        heavy_atoms_list.append(heavy_atoms)

    # Create dataframe and save as CSV
    df = pd.DataFrame({
        'MolecularID': mol_ids,
        'SMILES': smiles_list,
        'HeavyAtoms': heavy_atoms_list
    })
    df.to_csv(output_csv_path, index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python sdf2smileCSV.py input.sdf output.csv")
    else:
        input_sdf = sys.argv[1]
        output_csv = sys.argv[2]
        sdf_to_csv(input_sdf, output_csv)
        print(f"CSV written to {output_csv}")

