import re
import sys
import numpy as np
import pandas as pd
import os

def parse_xvg(fname):
    data = []
    header_labels = []
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('@'):
                match = re.search(r'(x|y)axis label\s+"(.+)"', line)
                if match:
                    header_labels.append(match.group(2).strip())
            elif line.startswith('#'):
                continue
            else:
                if line.strip() == '':
                    continue
                data.append(line.strip().split())
    arr = np.array(data, dtype=float)
    return arr, header_labels

def xvg_to_csv(xvg_file, csv_file=None, digits=7):
    arr, headers = parse_xvg(xvg_file)
    n_cols = arr.shape[1]
    if not headers or len(headers) != n_cols:
        headers = [f'Col{i+1}' for i in range(n_cols)]
    df = pd.DataFrame(arr, columns=headers)
    if not csv_file:
        base = os.path.splitext(xvg_file)[0]
        csv_file = base + '.csv'
    df.to_csv(csv_file, float_format=f'%.{digits}f', index=False)
    print(f'Converted {xvg_file} to {csv_file} with {n_cols} columns and {digits} digits.')

def batch_xvg_to_csv(folder_path, digits=7):
    for filename in os.listdir(folder_path):
        if filename.endswith('.xvg'):
            xvg_path = os.path.join(folder_path, filename)
            csv_path = os.path.splitext(xvg_path)[0] + '.csv'
            xvg_to_csv(xvg_path, csv_path, digits)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Convert all .xvg files in a folder to .csv format.")
    parser.add_argument('folder', type=str, help='Input folder containing XVG files')
    parser.add_argument('--digits', type=int, default=7, help='Decimal digits in output')
    args = parser.parse_args()
    batch_xvg_to_csv(args.folder, args.digits)

