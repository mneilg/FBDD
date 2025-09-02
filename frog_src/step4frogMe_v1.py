import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors, rdMolDescriptors
from tqdm import tqdm

# Check for correct usage
if len(sys.argv) != 3:
    print("Usage: python script.py input.tsv output.tsv")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

def rule_of_3(mol):
    mw = Descriptors.MolWt(mol)
    logP = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    return (150 <= mw <= 300 and logP <= 3 and hbd <= 3 and hba <= 3 and rot_bonds <= 3)

def generate_fragments(mol):
    try:
        brics_frags = BRICS.BRICSDecompose(mol)
        valid_frags = set()
        for frag_smiles in brics_frags:
            frag_mol = Chem.MolFromSmiles(frag_smiles)
            if frag_mol and rule_of_3(frag_mol):
                valid_frags.add(Chem.MolToSmiles(frag_mol))
        return valid_frags
    except Exception:
        return set()

df = pd.read_csv(input_file, sep='\t')
fragment_set = set()

for idx, row in tqdm(df.iterrows(), total=len(df)):
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        fragments = generate_fragments(mol)
        fragment_set.update(fragments)

frag_df = pd.DataFrame({'Fragment_SMILES': list(fragment_set)})
frag_df.to_csv(output_file, sep='\t', index=False)

print(f"Generated {len(fragment_set)} unique Rule-of-3 fragments.")


