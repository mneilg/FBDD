import sys
from rdkit import Chem
import matplotlib.pyplot as plt

def extract_property_list(sdf_file, property_name):
    values = []
    for mol in Chem.SDMolSupplier(sdf_file):
        if mol is None:
            continue
        try:
            val = float(mol.GetProp(property_name))
            values.append(val)
        except KeyError:
            continue  # Skip if property is missing
    return values

def plot_histogram(values, xlabel, ylabel, title, out_png):
    plt.figure(figsize=(8, 6))
    plt.hist(values, bins=30, color='skyblue', edgecolor='black')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_props.py input.sdf")
        sys.exit(1)
    sdf_file = sys.argv[1]

    # Extract and plot Molecular Weight
    mw_values = extract_property_list(sdf_file, 'MolecularWeight')
    plot_histogram(
        mw_values,
        xlabel='Molecular Weight (g/mol)',
        ylabel='Compound Count',
        title='Distribution of Molecular Weights',
        out_png='molecular_weight_histogram.png'
    )

    # Extract and plot logP
    logp_values = extract_property_list(sdf_file, 'LogP')
    plot_histogram(
        logp_values,
        xlabel='logP',
        ylabel='Compound Count',
        title='Distribution of logP Values',
        out_png='logp_histogram.png'
    )

