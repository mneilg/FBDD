import sys
from rdkit import Chem

def main():
    if len(sys.argv) != 3:
        print("Usage: python addSMILES.py input.sdf output.sdf")
        sys.exit(1)

    input_sdf = sys.argv[1]
    output_sdf = sys.argv[2]

    supplier = Chem.SDMolSupplier(input_sdf)
    writer = Chem.SDWriter(output_sdf)

    for mol in supplier:
        if mol is not None:
            # Convert to SMILES
            smiles = Chem.MolToSmiles(mol)
            # Add SMILES as a property
            mol.SetProp("SMILES", smiles)
            # Write to output SDF (all other properties are preserved)
            writer.write(mol)

    writer.close()

if __name__ == "__main__":
    main()

