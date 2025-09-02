import sys
from rdkit import Chem

def main():
    if len(sys.argv) != 2:
        print("Usage: python view_headers.py input.sdf")
        sys.exit(1)

    input_sdf = sys.argv[1]
    supplier = Chem.SDMolSupplier(input_sdf)
    headers = set()

    for mol in supplier:
        if mol is not None:
            for prop in mol.GetPropNames():
                headers.add(prop)

    print("Column headers (property names) in your SDF file:")
    for header in sorted(headers):
        print(header)

if __name__ == "__main__":
    main()

