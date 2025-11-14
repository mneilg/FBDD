import sys
import os
import argparse
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from tqdm import tqdm
import pandas as pd

def functional_group_counts(mol):
    fg_smarts = {
        'Amine': '[NX3;H2,H1;!$(NC=O)]',              # primary/secondary amines
        'Imine': '[CX3]=[NX2]',                       # C=N
        'Aldehyde': '[CX3H1](=O)[#6]',                # R-CHO
        'Ketone': '[CX3](=O)[#6]',                    # R2C=O
        'Ether': '[OD2]([#6])[#6]',                   # R-O-R
        'Carboxylic_Acid': 'C(=O)[OX2H1]',            # COOH
        'Sulfide': '[#16X2H0][#6]',                   # R-S-R
        'Thiol': '[#16X2H1]',                         # R-SH
    }
    return {name: len(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)))
            for name, smarts in fg_smarts.items()}

def calculate_ring_properties(mol):
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    has_complex_ring = any(len(ring) > 6 for ring in ring_sizes)
    return {
        'Ring_Count': len(ring_sizes),
        'Ring_Sizes': ','.join(map(str, ring_sizes)) if ring_sizes else '',
        'Complex_Rings': int(has_complex_ring)
    }

def element_counts(mol):
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    return {
        'Fluorine': symbols.count('F'),
        'Chlorine': symbols.count('Cl'),
        'Bromine': symbols.count('Br'),
        'Iodine': symbols.count('I'),
        'Sulfur': symbols.count('S')
    }

def hybridization_counts(mol):
    sp3 = sum(1 for atom in mol.GetAtoms() if atom.GetHybridization().name == 'SP3' and atom.GetSymbol() == 'C')
    sp2 = sum(1 for atom in mol.GetAtoms() if atom.GetHybridization().name == 'SP2' and atom.GetSymbol() == 'C')
    return {'sp3_Carbons': sp3, 'sp2_Carbons': sp2}

def calculate_descriptors(mol):
    desc = {}
    desc['TPSA'] = rdMolDescriptors.CalcTPSA(mol)
    desc['LogP'] = Descriptors.MolLogP(mol)
    desc['Rotatable_Bonds'] = Descriptors.NumRotatableBonds(mol)
    desc['Aromatic_Groups'] = Descriptors.NumAromaticRings(mol)

    desc.update(hybridization_counts(mol))
    desc.update(element_counts(mol))
    desc.update(calculate_ring_properties(mol))
    desc.update(functional_group_counts(mol))
    
    return desc

def main():
    parser = argparse.ArgumentParser(description="Extract physicochemical properties from an SDF file to a CSV file.")
    parser.add_argument("input_sdf", help="Path to input SDF file")
    parser.add_argument("output_csv", help="Path to output CSV file")
    args = parser.parse_args()

    if not os.path.isfile(args.input_sdf):
        print(f"Error: Input file '{args.input_sdf}' does not exist.")
        sys.exit(1)

    try:
        suppl = Chem.SDMolSupplier(args.input_sdf)
        if suppl is None or len(suppl) == 0:
            raise ValueError("No valid molecules found in SDF file.")
    except Exception as e:
        print(f"Error reading SDF file: {e}")
        sys.exit(2)

    data = []

    for idx, mol in enumerate(tqdm(suppl, desc="Processing molecules")):
        if mol is None:
            print(f"Warning: Molecule {idx + 1} is invalid and will be skipped.")
            continue

        entry = {
            'Molecule_ID': mol.GetProp("_Name") if mol.HasProp("_Name") else f"MOL_{idx+1}",
            'SMILES': Chem.MolToSmiles(mol)
        }
        entry.update(calculate_descriptors(mol))
        data.append(entry)

    if not data:
        print("Error: No valid molecules processed.")
        sys.exit(3)

    try:
        df = pd.DataFrame(data)
        df.to_csv(args.output_csv, index=False)
        print(f"CSV successfully written to {args.output_csv}")
    except Exception as e:
        print(f"Error writing CSV file: {e}")
        sys.exit(4)

if __name__ == "__main__":
    main()

