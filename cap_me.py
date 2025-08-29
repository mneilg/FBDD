import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolops

def cap_one_dummy(mol, dummy_idx):
    atom = mol.GetAtomWithIdx(dummy_idx)
    neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
    # Check if dummy is in a ring
    is_ring = atom.IsInRing()
    em = Chem.EditableMol(mol)
    em.RemoveAtom(dummy_idx)
    mol2 = em.GetMol()
    # After removal, neighbor indices may have shifted.
    # Use map to relate old to new indices
    map_idx = {}
    cnt = 0
    for i in range(mol.GetNumAtoms()):
        if i != dummy_idx:
            map_idx[i] = cnt
            cnt += 1
    # Attach appropriate cap (H or CH3) to each neighbor
    rw = Chem.RWMol(mol2)
    for nb_idx in neighbors:
        real_nb_idx = map_idx[nb_idx]
        if is_ring:
            frag = Chem.MolFromSmiles('C')
        else:
            frag = Chem.MolFromSmiles('[H]')
        new_atom_idx = rw.AddAtom(frag.GetAtomWithIdx(0))
        rw.AddBond(real_nb_idx, new_atom_idx, Chem.rdchem.BondType.SINGLE)
    return rw.GetMol()

def robust_cap_all_dummies(smiles_in):
    # Returns sanitized, capped SMILES or None if failed
    try:
        mol = Chem.MolFromSmiles(smiles_in, sanitize=False)
        if mol is None:
            return None
        # Update ring info, in case it's needed
        mol.UpdatePropertyCache(strict=False)
        rdmolops.GetSymmSSSR(mol)
        dummy_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == '*']
        # Cap all dummies, one at a time
        while dummy_idxs:
            # Always cap the lowest index first
            dummy_idxs = sorted(dummy_idxs)
            mol = cap_one_dummy(mol, dummy_idxs[0])
            # Update ring info after edit
            mol.UpdatePropertyCache(strict=False)
            rdmolops.GetSymmSSSR(mol)
            dummy_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == '*']
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception as e:
        return None

def cap_fragments(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t")
    rows = []
    for smi in df['Fragment_SMILES']:
        capped = robust_cap_all_dummies(str(smi))
        if capped is not None:
            rows.append([capped])
    outdf = pd.DataFrame(rows, columns=['Fragment_SMILES'])
    outdf.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python cap_me.py input.tsv output.tsv")
        sys.exit(1)
    cap_fragments(sys.argv[1], sys.argv[2])

