import sys
from rdkit import Chem
from rdkit.Chem import PandasTools, inchi
from rdkit.Chem import SDWriter
from molvs import Standardizer

if len(sys.argv) != 3:
    print("Usage: python process_sdf.py <input.sdf> <output.sdf>")
    sys.exit(1)

input_sdf = sys.argv[1]
output_sdf = sys.argv[2]

# 1. Load SDF and add SMILES
frame = PandasTools.LoadSDF(input_sdf, smilesName='SMILES', molColName='Molecule', includeFingerprints=False)

# 2. Standardize (remove salts, canonicalize tautomers, validate)
standardizer = Standardizer()
frame['Molecule'] = frame['Molecule'].apply(lambda m: standardizer.standardize(m) if m else None)

# 3. Remove duplicates (InChIKey)
frame['InChIKey'] = frame['Molecule'].apply(lambda m: inchi.MolToInchiKey(m) if m else None)
frame = frame.drop_duplicates(subset='InChIKey')

# 4. Filter out invalid molecules (optional, but Standardizer already does this)
frame = frame[frame['Molecule'].notna()]

# 5. Save to SDF
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

