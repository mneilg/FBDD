import sys
from rdkit import Chem
from rdkit.Chem import PandasTools, inchi
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import rdMolStandardize
from rdkit.Chem import SDWriter

if len(sys.argv) != 3:
    print("Usage: python process_sdf.py <input.sdf> <output.sdf>")
    sys.exit(1)

input_sdf = sys.argv[1]
output_sdf = sys.argv[2]

# 1. Load SDF and add SMILES
frame = PandasTools.LoadSDF(input_sdf, smilesName='SMILES', molColName='Molecule', includeFingerprints=False)

# 2. Remove salts/counterions
remover = SaltRemover()
frame['Molecule'] = frame['Molecule'].apply(lambda m: remover.StripMol(m, dontRemoveEverything=True) if m else None)

# 3. Canonicalize tautomers
enumerator = rdMolStandardize.TautomerEnumerator()
frame['Molecule'] = frame['Molecule'].apply(lambda m: enumerator.Canonicalize(m) if m else None)

# 4. Remove duplicates (InChIKey)
frame['InChIKey'] = frame['Molecule'].apply(lambda m: inchi.MolToInchiKey(m) if m else None)
frame = frame.drop_duplicates(subset='InChIKey')

# 5. Validate structures
def is_valid_mol(mol):
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
        return True
    except:
        return False
frame = frame[frame['Molecule'].apply(is_valid_mol)]

# 6. Save to SDF
writer = SDWriter(output_sdf)
for idx, row in frame.iterrows():
    mol = row['Molecule']
    if mol is None:
        continue
    mol.SetProp('SMILES', row['SMILES'])
    mol.SetProp('InChIKey', row['InChIKey'])
    writer.write(mol)
writer.close()

print(f"Successfully processed and saved to {output_sdf}")

