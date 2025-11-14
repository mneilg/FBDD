import sys
import os
import argparse
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, MACCSkeys
from rdkit.Chem import AllChem
from tqdm import tqdm
import pandas as pd

def functional_group_counts(mol):
    fg_smarts = {
        'Alcohol': '[CX4][OH]',
        'Amine': '[NX3;H2,H1;!$(NC=O)]',
        'Carboxyl': '[CX3](=O)[OX2H1]',
        'Carbonyl': '[CX3]=O',
        'Imine': '[CX3]=[NX2]', # C=N
        'Thiol': '[SX2H]', # S-H
        'Sulfide': '[SX2][CX4]', # S-C
        'Halide': '[#6][F,Cl,Br,I]', # Alkyl halide
        'Amide': '[NX3][CX3](=O)[#6]', # N-C(=O)-C
    }
    fg_counts = {name: len(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)))
                 for name, smarts in fg_smarts.items()}
    return fg_counts

def conjugated_bonds_count(mol):
    # Counts bonds where both atoms are in conjugated systems
    return sum(1 for bond in mol.GetBonds() if bond.GetIsConjugated())

def calculate_descriptors(mol):
    desc = {}
    desc['MolWt'] = Descriptors.MolWt(mol)  # Molecular weight [g/mol]
    desc['LogP'] = Descriptors.MolLogP(mol) # Octanol-water logP
    desc['TPSA'] = rdMolDescriptors.CalcTPSA(mol)
    desc['Num_HBD'] = Descriptors.NumHDonors(mol)
    desc['Num_HBA'] = Descriptors.NumHAcceptors(mol)
    desc['Rotatable_Bonds'] = Descriptors.NumRotatableBonds(mol)
    desc['Conjugated_Bonds_Count'] = conjugated_bonds_count(mol)
    desc['Ring_Count'] = rdMolDescriptors.CalcNumRings(mol)
    desc['Aromatic_Rings'] = Descriptors.NumAromaticRings(mol)
    desc['Halogen_Count'] = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() in ['F', 'Cl', 'Br', 'I'])
    desc['Heavy_Atoms'] = Descriptors.HeavyAtomCount(mol)
    #desc['WienerIndex'] = rdMolDescriptors.CalcWienerIndex(mol)
    #desc['HosoyaIndex'] = rdMolDescriptors.CalcHosoyaIndex(mol)
    desc.update(functional_group_counts(mol))
    return desc

def main():
    parser = argparse.ArgumentParser(description="Calculate molecular descriptors from SDF file.")
    parser.add_argument("input_sdf", help="Path to input SDF file")
    parser.add_argument("output_csv", help="Path to output CSV file")
    args = parser.parse_args()

    if not os.path.isfile(args.input_sdf):
        print(f"Error: Input file '{args.input_sdf}' does not exist.")
        sys.exit(1)

    try:
        suppl = Chem.SDMolSupplier(args.input_sdf)
        if suppl is None or len(suppl) == 0:
            raise ValueError("No molecules found or input file could not be read as SDF.")
    except Exception as e:
        print(f"Error reading SDF file: {e}")
        sys.exit(2)

    data = []
    for idx, mol in enumerate(tqdm(suppl, desc="Processing molecules")):
    	if mol is None:
        	print(f"Warning: Molecule {idx+1} in SDF is invalid and will be skipped.")
        	continue

    	entry = {
        	'MolecularID': mol.GetProp("MolecularID") if mol.HasProp("MolecularID") else "",
        	'LogP': mol.GetProp("LogP") if mol.HasProp("LogP") else "",
        	'MolecularWeight': mol.GetProp("MolecularWeight") if mol.HasProp("MolecularWeight") else "",
        	'SMILES': Chem.MolToSmiles(mol)
    	}
    	entry.update(calculate_descriptors(mol))
    	data.append(entry)

    if not data:
        print("Error: No valid molecules processed. No output generated.")
        sys.exit(3)

    try:
        df = pd.DataFrame(data)
        df.to_csv(args.output_csv, index=False)
        print(f"CSV file successfully written to: {args.output_csv}")
    except Exception as e:
        print(f"Error writing CSV file: {e}")
        sys.exit(4)

if __name__ == "__main__":
    main()

