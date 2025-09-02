import sys
import pandas as pd
from rdkit import Chem

def main(infile, outfile):
    df = pd.read_csv(infile, sep='\t')
    if 'ID' not in df.columns:
        df['ID'] = [f'F{str(i+1).zfill(4)}' for i in range(len(df))]
    valid_lines = []
    for idx, row in df.iterrows():
        smi = str(row['Fragment_SMILES']).strip()
        id_ = str(row['ID']).strip()
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        has_unassigned = any(c[1] == '?' for c in chiral_centers)
        if has_unassigned:
            continue
        valid_lines.append(f"{smi} {id_}\n")
    with open(outfile, 'w') as f:
        f.writelines(valid_lines)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python cleanfrogs.py input.tsv output.smi")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

