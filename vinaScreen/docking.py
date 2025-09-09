#!/usr/bin/env python3
"""
docking.py - Run AutoDock Vina for a receptor-ligand pair using your config.txt
- Handles output file naming (concatenated receptor/ligand, numbered for duplicates)
- Parses affinity and RMSD from Vina output
- Returns status and error info for reporting
"""

import os
import subprocess
import re
import sys
from pathlib import Path

def generate_output_filename(receptor, ligand, output_dir="output"):
    rec_name = Path(receptor).stem
    lig_name = Path(ligand).stem
    base = f"{rec_name}_{lig_name}"
    filename = f"{base}.pdbqt"
    fullpath = os.path.join(output_dir, filename)
    # Add a numeric suffix if duplicate exists
    i = 1
    while os.path.exists(fullpath):
        filename = f"{base}_{i}.pdbqt"
        fullpath = os.path.join(output_dir, filename)
        i += 1
    return fullpath

def run_vina(receptor, ligand, config="config.txt", vina_exec="vina_1.2.7_linux_x86_64", output_dir="output"):
    if not os.path.exists(vina_exec):
        return {"success": False, "error": "vina_1.2.7_linux_x86_64 executable is not found in the current directory"}
    if not os.path.exists(receptor):
        return {"success": False, "error": f"Receptor file (`{receptor}`) not found"}
    if not os.path.exists(ligand):
        return {"success": False, "error": f"Ligand file (`{ligand}`) not found"}
    if not os.path.exists(config):
        return {"success": False, "error": "config.txt file is not found in the current directory"}
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_file = generate_output_filename(receptor, ligand, output_dir)
    cmd = [
        f"./{vina_exec}",
        "--receptor", receptor,
        "--ligand", ligand,
        "--config", config,
        "--out", output_file
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        out = result.stdout + result.stderr
        if result.returncode != 0:
            return {"success": False, "error": out}

        # Parse affinity and RMSD from output
        affinity, rmsd_lb = None, None
        # Standard Vina output: lines like "   1       -8.5      0.000      0.000"
        for line in out.split('\n'):
            if re.match(r'^\s*1\s+', line):
                parts = line.split()
                if len(parts) >= 4:
                    affinity = float(parts[1])
                    rmsd_lb = float(parts[2])
        # If not found, report parsing error
        if affinity is None:
            return {"success": False, "error": "Could not parse affinity from Vina output"}
        return {
            "success": True,
            "output_file": output_file,
            "affinity": affinity,
            "rmsd_lb": rmsd_lb
        }
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "Vina execution timed out"}
    except Exception as exc:
        return {"success": False, "error": f"Unexpected error: {exc}"}

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python docking.py <receptor_file> <ligand_file>")
        print("Example: python docking.py receptor/protein.pdbqt ligand/compound.pdbqt")
        sys.exit(1)
    receptor, ligand = sys.argv[1], sys.argv[2]
    result = run_vina(receptor, ligand)
    if result["success"]:
        print(f"✓ Docking successful!\n  Output: {result['output_file']}\n  Affinity: {result['affinity']} kcal/mol\n  RMSD_lb: {result['rmsd_lb']}")
        sys.exit(0)
    else:
        print(f"✗ Docking failed: {result['error']}")
        sys.exit(1)

