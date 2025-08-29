import sys
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops, SDWriter

def cap_one_dummy(mol, dummy_idx):
    atom = mol.GetAtomWithIdx(dummy_idx)
    neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
    is_ring = atom.IsInRing()

    # Mark all atom indices before removal for mapping
    num_atoms = mol.GetNumAtoms()
    idx_map = {i: i for i in range(num_atoms)}  # original to current index

    em = Chem.EditableMol(mol)
    em.RemoveAtom(dummy_idx)
    mol2 = em.GetMol()

    # After atom removal, indices shift.
    # We need to build a mapping from old indices to new, skipping the dummy atom
    map_idx = {}
    cnt = 0
    for i in range(num_atoms):
        if i != dummy_idx:
            map_idx[i] = cnt
            cnt += 1
    # Only keep neighbor indices that exist post-removal
    valid_neighbors = [nb for nb in neighbors if nb in map_idx]

    rw = Chem.RWMol(mol2)
    for nb_idx in valid_neighbors:
        real_nb_idx = map_idx[nb_idx]
        if is_ring:
            frag = Chem.MolFromSmiles('C')
        else:
            frag = Chem.MolFromSmiles('[H]')
        new_atom_idx = rw.AddAtom(frag.GetAtomWithIdx(0))
        rw.AddBond(real_nb_idx, new_atom_idx, Chem.rdchem.BondType.SINGLE)
    cap_mol = rw.GetMol()

    # Sanitize, remove unconnected Hs if present
    try:
        Chem.SanitizeMol(cap_mol)
    except Exception:
        # Try to remove fragments that are lone H and sanitize again
        fragments = Chem.GetMolFrags(cap_mol, asMols=True, sanitizeFrags=False)
        # Remove any fragment that's only H
        fragments = [frag for frag in fragments if frag.GetNumAtoms() > 1 or frag.GetAtomWithIdx(0).GetSymbol() != 'H']
        if fragments:
            cap_mol = fragments  # take main frag
            Chem.SanitizeMol(cap_mol)
        else:
            return None
    return cap_mol

def robust_cap_all_dummies(mol):
    try:
        mol.UpdatePropertyCache(strict=False)
        rdmolops.GetSymmSSSR(mol)
        dummy_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() in ['*', 'R']]
        while dummy_idxs:
            # Cap the dummy at the lowest-index each time: REMAP indices each loop!
            dummy_idxs = sorted(dummy_idxs)
            mol = cap_one_dummy(mol, dummy_idxs)
            if mol is None:
                return None
            mol.UpdatePropertyCache(strict=False)
            rdmolops.GetSymmSSSR(mol)
            dummy_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() in ['*', 'R']]
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None

def gen_3d(mol):
    try:
        mol_3d = Chem.AddHs(mol)
        res = AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
        if res == 0:
            AllChem.UFFOptimizeMolecule(mol_3d)
        Chem.SanitizeMol(mol_3d)
        return mol_3d
    except Exception:
        return None

def process_sdf(input_sdf, output_sdf, id_prefix="Frag"):
    suppl = Chem.SDMolSupplier(input_sdf, removeHs=False)
    writer = SDWriter(output_sdf)
    count = 1
    for mol in suppl:
        if mol is None:
            continue
        capped_mol = robust_cap_all_dummies(mol)
        if capped_mol is None:
            continue
        mol_3d = gen_3d(capped_mol)
        if mol_3d is None:
            continue
        mol_id = f"{id_prefix}{count:05d}"
        mol_3d.SetProp("MolecularID", mol_id)
        mol_3d.SetProp("CappedSMILES", Chem.MolToSmiles(mol_3d, isomericSmiles=True))
        writer.write(mol_3d)
        count += 1
    writer.close()
    print(f"Saved {count-1} capped, 3D molecules to {output_sdf}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python cap_and_embed_sdf.py input.sdf output.sdf")
        sys.exit(1)
    process_sdf(sys.argv[1], sys.argv[2])

