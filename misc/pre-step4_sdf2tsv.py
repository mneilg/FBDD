import sys
from rdkit import Chem
from rdkit.Chem import Descriptors

def extract_heavy_molecules(input_sdf, output_tsv):
    supplier = Chem.SDMolSupplier(input_sdf)
    with open(output_tsv, 'w') as f:
        # Write header
        f.write("InChIKey\tLogP\tMolecularID\tMolecularWeight\tSMILES\n")
        for mol in supplier:
            if mol is None:
                continue
            try:
                mol_weight = float(mol.GetProp("MolecularWeight"))
            except KeyError:
                mol_weight = Descriptors.MolWt(mol)
            # Size filter removed: all molecules are included
            try:
                inchikey = Chem.MolToInchiKey(mol)
            except Exception:
                inchikey = ""
            try:
                logp = mol.GetProp("LogP")
            except KeyError:
                logp = ""
            try:
                mol_id = mol.GetProp("MolecularID")
            except KeyError:
                mol_id = ""
            try:
                smiles = Chem.MolToSmiles(mol)
            except Exception:
                smiles = ""
            f.write(f"{inchikey}\t{logp}\t{mol_id}\t{mol_weight:.2f}\t{smiles}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_heavy.py input.sdf output.tsv")
        sys.exit(1)
    input_sdf = sys.argv[1]
    output_tsv = sys.argv[2]
    extract_heavy_molecules(input_sdf, output_tsv)

