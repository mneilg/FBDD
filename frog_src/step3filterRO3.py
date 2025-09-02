import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def is_ro3(mol):
    """Checks if molecule satisfies Rule of Three (RO3)."""
    try:
        mw = float(mol.GetProp('MolecularWeight'))
    except Exception:
        mw = Descriptors.MolWt(mol)
    try:
        logp = float(mol.GetProp('LogP'))
    except Exception:
        logp = Crippen.MolLogP(mol)
    hba = Chem.rdMolDescriptors.CalcNumHBA(mol)
    hbd = Chem.rdMolDescriptors.CalcNumHBD(mol)
    rot_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
    return (
        mw <= 300 and
        logp <= 3 and
        hba <= 3 and
        hbd <= 3 and
        rot_bonds <= 3
    )

def filter_ro3(input_sdf, output_pass, output_fail):
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    writer_pass = Chem.SDWriter(output_pass)
    writer_fail = Chem.SDWriter(output_fail)
    count_pass, count_fail, count_problem = 0, 0, 0

    for mol in supplier:
        if mol is None:
            count_problem += 1
            continue
        try:
            # Ensure molecule is sanitized and valid
            Chem.SanitizeMol(mol)
            if is_ro3(mol):
                writer_pass.write(mol)
                count_pass += 1
            else:
                writer_fail.write(mol)
                count_fail += 1
        except Exception:
            count_problem += 1
            continue
    writer_pass.close()
    writer_fail.close()
    print(f"Processed: {count_pass} passed, {count_fail} failed, {count_problem} problematic molecules removed.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_ro3.py input.sdf passed.sdf failed.sdf")
        sys.exit(1)
    input_sdf = sys.argv[1]
    output_pass = sys.argv[2]
    output_fail = sys.argv[3]
    filter_ro3(input_sdf, output_pass, output_fail)

