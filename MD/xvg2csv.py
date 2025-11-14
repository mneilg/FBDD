import re
import sys
import numpy as np
import pandas as pd
import os

def parse_xvg(fname):
    # Read file, skip comment/header lines (@, #)
    data = []
    header_labels = []
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('@'):
                # Attempt to get labels from xaxis/yaxis legends
                match = re.search(r'(x|y)axis label\s+"(.+)"', line)
                if match:
                    header_labels.append(match.group(2).strip())
            elif line.startswith('#'):
                continue
            else:
                # Read data line
                if line.strip() == '':
                    continue
                data.append(line.strip().split())
    arr = np.array(data, dtype=float)
    return arr, header_labels

def xvg_to_csv(xvg_file, csv_file=None, digits=7):
    arr, headers = parse_xvg(xvg_file)
    # Generate generic headers if missing
    n_cols = arr.shape[1]
    if not headers or len(headers) != n_cols:
        headers = [f'Col{i+1}' for i in range(n_cols)]
    df = pd.DataFrame(arr, columns=headers)
    # Output file logic
    if not csv_file:
        base = os.path.splitext(xvg_file)[0]
        csv_file = base + '.csv'
    df.to_csv(csv_file, float_format=f'%.{digits}f', index=False)
    print(f'Converted {xvg_file} to {csv_file} with {n_cols} columns and {digits} digits.')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Convert .xvg file to .csv, handling multi-column data.")
    parser.add_argument('xvg', type=str, help='Input XVG file')
    parser.add_argument('--csv', type=str, default=None, help='Output CSV file (optional)')
    parser.add_argument('--digits', type=int, default=7, help='Decimal digits in output')
    args = parser.parse_args()
    xvg_to_csv(args.xvg, args.csv, args.digits)

