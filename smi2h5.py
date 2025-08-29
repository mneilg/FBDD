import sys
import h5py
import numpy as np

if len(sys.argv) != 3:
    print("Usage: python smi2h5.py input.smi output.h5")
    sys.exit(1)

input_smi = sys.argv[1]
output_h5 = sys.argv[2]

# Read SMILES from the input .smi file
with open(input_smi, 'r') as f:
    smiles = [line.strip().split()[0] for line in f if line.strip()]

# Convert to numpy array of ASCII-encoded bytes
smiles_bytes = np.array([s.encode('ascii') for s in smiles])

# Write to HDF5 with dataset name 'frag_smiles'
with h5py.File(output_h5, 'w') as f:
    f.create_dataset('frag_smiles', data=smiles_bytes)

