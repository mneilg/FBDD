import sys
import os
import tempfile
from rdkit import Chem
from rdkit.Chem import SDWriter
from tqdm import tqdm

def get_smiles(mol):
    try:
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return ""

def validate_with_openbabel(mol):
    try:
        block = Chem.MolToMolBlock(mol)
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.mol') as temp_mol:
            temp_mol.write(block)
            temp_mol.flush()
            temp_mol_name = temp_mol.name
        with tempfile.NamedTemporaryFile(mode='r', delete=False, suffix='.smi') as temp_smi:
            temp_smi_name = temp_smi.name
        cmd = f'obabel -imol "{temp_mol_name}" -osmi -O "{temp_smi_name}"'
        ret = os.system(cmd)
        ok = (ret == 0) and os.path.exists(temp_smi_name) and os.path.getsize(temp_smi_name) > 0
        os.remove(temp_mol_name)
        os.remove(temp_smi_name)
        return ok
    except Exception:
        return False

def read_valid_mols_sdf(sdf_file, id_prefix, start_id):
    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
    suppl = [mol for mol in suppl if mol is not None]
    count = start_id
    written = 0
    with tqdm(total=len(suppl), desc=f"Processing {os.path.basename(sdf_file)}") as pbar:
        for mol in suppl:
            pbar.update(1)
            try:
                Chem.SanitizeMol(mol)
                frags = Chem.GetMolFrags(mol)
                if len(frags) > 1:
                    largest = max(frags, key=len)
                    mol = Chem.PathToSubmol(mol, largest)
                    Chem.SanitizeMol(mol)
                mol.SetProp("MolecularID", f"{id_prefix}{count:06d}")
                smiles = get_smiles(mol)
                mol.SetProp("SMILES", smiles)
                if not validate_with_openbabel(mol):
                    continue
                count += 1
                written += 1
                yield mol
            except Exception:
                continue
   # print(f"{sdf_file}: Done. {written} valid molecules written.")

def combine_to_sdf(sdf_file1, sdf_file2, output_file):
    writer = SDWriter(output_file)
    count = 1
    for mol in read_valid_mols_sdf(sdf_file1, "MOL_", count):
        writer.write(mol)
        count += 1
    for mol in read_valid_mols_sdf(sdf_file2, "MOL_", count):
        writer.write(mol)
        count += 1
    writer.close()
    print(f"Combined SDF generated: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python combine_Ro3Frogs.py input1.sdf input2.sdf output.sdf")
        sys.exit(1)
    combine_to_sdf(sys.argv[1], sys.argv[2], sys.argv[3])

