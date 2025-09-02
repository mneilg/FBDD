from rdkit import Chem
import os

def rdkit_merge_with_id(input_dir, output_file):
    writer = Chem.SDWriter(output_file)
    for filename in os.listdir(input_dir):
        if filename.endswith('.mol'):
            filepath = os.path.join(input_dir, filename)
            mol = Chem.MolFromMolFile(filepath)
            if mol:
                # Use filename (without extension) as molecule ID
                mol.SetProp("ID", os.path.splitext(filename)[0])
                writer.write(mol)
    writer.close()

if __name__ == "__main__":
    import sys
    rdkit_merge_with_id(sys.argv[1], sys.argv[2])

