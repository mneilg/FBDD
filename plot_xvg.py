import matplotlib.pyplot as plt
import re
import argparse

def plot_xvg(filename, title, save_file):
    x = []
    y = []
    xlabel = 'X'
    ylabel = 'Y'
    # Read file, extract data and axis labels
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('@'):
                match_x = re.match(r'@ *xaxis *label *\"(.+)\"', line)
                match_y = re.match(r'@ *yaxis *label *\"(.+)\"', line)
                if match_x:
                    xlabel = match_x.group(1)
                if match_y:
                    ylabel = match_y.group(1)
                continue
            if line.startswith('#'):
                continue
            cols = line.strip().split()
            if len(cols) >= 2:
                x.append(float(cols[0]))
                y.append(float(cols[1]))
    plt.figure(figsize=(8, 5))
    plt.plot(x, y, label=filename)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_file, dpi=300, bbox_inches='tight')
    print(f"PNG image saved as {save_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot GROMACS XVG file and save as PNG")
    parser.add_argument('-f', '--file', required=True, help='Input XVG filename')
    parser.add_argument('-t', '--title', required=True, help='Plot title')
    parser.add_argument('-s', '--save', required=True, help='Output PNG filename')
    args = parser.parse_args()
    plot_xvg(args.file, args.title, args.save)

