import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from tqdm import tqdm
import networkx as nx


def mol_to_nx(rdkit_mol):
    G = nx.Graph()
    for bond in rdkit_mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return G

def expand_fingerprint(fp, length):
    return [int(fp.GetBit(i)) for i in range(length)]

def main(input_sdf, output_csv):
    suppl = Chem.SDMolSupplier(input_sdf)
    data = []
    for mol in tqdm(suppl, desc="Processing molecules"):
        if mol is None:
            continue
        entry = {
            'MolecularID': mol.GetProp("MolecularID") if mol.HasProp("MolecularID") else "",
            'SMILES': Chem.MolToSmiles(mol)
        }
        ecfp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        ecfp_list = expand_fingerprint(ecfp, 2048)
        for i, bit in enumerate(ecfp_list):
            entry[f"ECFP4_{i}"] = bit
        maccs = MACCSkeys.GenMACCSKeys(mol)
        maccs_list = list(maccs)
        for i, bit in enumerate(maccs_list):
            entry[f"MACCS_{i}"] = bit

        nx_graph = mol_to_nx(mol)
        try:
            entry["WienerIndex"] = wiener_index(nx_graph)
        except Exception:
            entry["WienerIndex"] = ""
        # Hosoya index is not supported out of the box in modern Python libraries

        data.append(entry)
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    print(f"Saved output to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate ECFP, MACCS, and Wiener indices from SDF (NetworkX)")
    parser.add_argument("input_sdf", help="Input SDF file path")
    parser.add_argument("output_csv", help="Output CSV file path")
    args = parser.parse_args()
    main(args.input_sdf, args.output_csv)

