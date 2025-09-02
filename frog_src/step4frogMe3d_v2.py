import sys
from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors, rdMolDescriptors, AllChem, SDWriter
from tqdm import tqdm
import pandas as pd

if len(sys.argv) != 3:
    print("Usage: python step4frogMe_3D.py input.sdf output_prefix")
    sys.exit(1)

input_sdf = sys.argv[1]
output_prefix = sys.argv[2]

def rule_of_3(mol):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    return (150 <= mw <= 300 and logp <= 3 and hbd <= 3 and hba <= 3 and rot_bonds <= 3)

suppl = Chem.SDMolSupplier(input_sdf, removeHs=False)
all_fragments_smiles = set()
all_fragments_mols = []

for mol in tqdm(suppl):
    if mol is None:
        continue
    try:
        brics_frags = BRICS.BRICSDecompose(mol)
    except Exception:
        brics_frags = []
    for frag_smiles in brics_frags:
        frag_mol = Chem.MolFromSmiles(frag_smiles)
        if frag_mol is None:
            continue
        if rule_of_3(frag_mol):
            # Only keep fragments with at least one dummy atom ('[*]')
            if any(atom.GetAtomicNum() == 0 for atom in frag_mol.GetAtoms()):
                frag_smiles_std = Chem.MolToSmiles(frag_mol)
                if frag_smiles_std not in all_fragments_smiles:
                    all_fragments_smiles.add(frag_smiles_std)
                    all_fragments_mols.append(frag_mol)



# Write fragments—3D SDF
writer = SDWriter(f"{output_prefix}_3dFrogs.sdf")
for m in all_fragments_mols:
    writer.write(m)
writer.close()

# Also export as SMILES for 2D processing
frag_smiles_list = [Chem.MolToSmiles(m) for m in all_fragments_mols]
frag_df = pd.DataFrame({"Fragment_SMILES": frag_smiles_list})
frag_df.to_csv(f"{output_prefix}_FroggySmiles.tsv", sep="\t", index=False)

print(f"Generated {len(all_fragments_mols)} unique Rule-of-3 fragments with 3D, saved as {output_prefix}_3DFrogs.sdf and {output_prefix}_FroggySmiles.tsv")

