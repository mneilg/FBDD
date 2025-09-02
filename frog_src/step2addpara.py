import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def enhance_sdf_properties(input_sdf, output_sdf):
    supplier = Chem.SDMolSupplier(input_sdf)
    writer = Chem.SDWriter(output_sdf)
    
    for idx, mol in enumerate(supplier):
        if mol is None:
            continue
        
        # Calculate molecular weight
        mol_weight = Descriptors.MolWt(mol)
        mol.SetProp('MolecularWeight', f"{mol_weight:.4f}")
        
        # Calculate Crippen logP
        logp = Crippen.MolLogP(mol)
        mol.SetProp('LogP', f"{logp:.4f}")
        
        # Generate sequential molecular ID
        mol_id = f"q{idx+1:05d}"
        mol.SetProp('MolecularID', mol_id)
        
        writer.write(mol)
    
    writer.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python addpara.py input.sdf output.sdf")
        sys.exit(1)
    input_sdf = sys.argv[1]
    output_sdf = sys.argv[2]
    enhance_sdf_properties(input_sdf, output_sdf)

