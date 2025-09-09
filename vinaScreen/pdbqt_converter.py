#!/usr/bin/env python3
"""
OpenBabel PDBQT Conversion Script for AutoDock Vina
Converts molecular libraries (.smi, .csv, .tsv, .sdf, .mol) to PDBQT format
Author: Chemistry Library Converter
Version: 1.0
"""

import subprocess
import sys
import os
import tempfile
import pandas as pd
import argparse
from pathlib import Path
import time

TIMEOUT_SECS = 20  # Default timeout for OpenBabel conversion

class PDBQTConverter:
    def __init__(self, timeout=TIMEOUT_SECS):
        self.timeout = timeout
        self.stats = {
            'total': 0,
            'converted': 0,
            'skipped': 0,
            'errors': 0
        }

    def detect_file_format(self, filepath):
        """Detect input file format based on extension"""
        extension = Path(filepath).suffix.lower()
        if extension in ['.smi', '.smiles']:
            return 'smiles'
        elif extension in ['.csv']:
            return 'csv'
        elif extension in ['.tsv']:
            return 'tsv'
        elif extension in ['.sdf']:
            return 'sdf'
        elif extension in ['.mol']:
            return 'mol'
        else:
            raise ValueError(f"Unsupported file format: {extension}")

    def read_smiles_file(self, filepath):
        """Read SMILES file and return list of (smiles, molecule_id) tuples"""
        molecules = []
        try:
            with open(filepath, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    parts = line.split()
                    if len(parts) >= 2:
                        smiles, mol_id = parts[0], parts[1]
                    elif len(parts) == 1:
                        smiles, mol_id = parts[0], f"mol_{line_num:04d}"
                    else:
                        continue

                    molecules.append({
                        'smiles': smiles,
                        'mol_id': mol_id,
                        'mol_block': None
                    })
        except Exception as e:
            print(f"Error reading SMILES file: {e}")
            return []

        return molecules

    def read_csv_tsv_file(self, filepath, delimiter=','):
        """Read CSV/TSV file with SMILES and/or MOL columns"""
        molecules = []
        try:
            df = pd.read_csv(filepath, delimiter=delimiter)

            # Try to identify relevant columns
            smiles_cols = [col for col in df.columns if 'smiles' in col.lower()]
            mol_cols = [col for col in df.columns if 'mol' in col.lower() and 'smiles' not in col.lower()]
            id_cols = [col for col in df.columns if any(x in col.lower() for x in ['id', 'name', 'title'])]

            smiles_col = smiles_cols[0] if smiles_cols else None
            mol_col = mol_cols[0] if mol_cols else None
            id_col = id_cols[0] if id_cols else None

            if not smiles_col and not mol_col:
                # Assume first column is SMILES if no clear column found
                smiles_col = df.columns[0]
                print(f"Warning: No clear SMILES/MOL column found, using '{smiles_col}'")

            for idx, row in df.iterrows():
                mol_id = row[id_col] if id_col else f"mol_{idx+1:04d}"
                smiles = row[smiles_col] if smiles_col and pd.notna(row[smiles_col]) else None
                mol_block = row[mol_col] if mol_col and pd.notna(row[mol_col]) else None

                if smiles or mol_block:
                    molecules.append({
                        'smiles': smiles,
                        'mol_id': str(mol_id),
                        'mol_block': mol_block
                    })

        except Exception as e:
            print(f"Error reading CSV/TSV file: {e}")
            return []

        return molecules

    def read_sdf_file(self, filepath):
        """Read SDF file and extract molecules"""
        molecules = []
        try:
            with open(filepath, 'r') as f:
                content = f.read()

            # Split by molecule delimiter
            mol_blocks = content.split('$$$$')

            for i, mol_block in enumerate(mol_blocks):
                mol_block = mol_block.strip()
                if not mol_block:
                    continue

                # Extract molecule name from first line or use default
                lines = mol_block.split('\n')
                mol_id = lines[0].strip() if lines[0].strip() else f"mol_{i+1:04d}"

                molecules.append({
                    'smiles': None,
                    'mol_id': mol_id,
                    'mol_block': mol_block
                })

        except Exception as e:
            print(f"Error reading SDF file: {e}")
            return []

        return molecules

    def read_mol_file(self, filepath):
        """Read single MOL file"""
        try:
            with open(filepath, 'r') as f:
                mol_block = f.read().strip()

            mol_id = Path(filepath).stem
            return [{
                'smiles': None,
                'mol_id': mol_id,
                'mol_block': mol_block
            }]

        except Exception as e:
            print(f"Error reading MOL file: {e}")
            return []

    def convert_molecule_to_pdbqt(self, molecule, output_path, tmpdir):
        """Convert a single molecule to PDBQT format"""
        mol_id = molecule['mol_id']
        smiles = molecule['smiles']
        mol_block = molecule['mol_block']

        # Prepare temporary files
        temp_input = os.path.join(tmpdir, f"{mol_id}_input")
        temp_pdbqt = os.path.join(tmpdir, f"{mol_id}_temp.pdbqt")

        # Try MOL block first (preferred), then SMILES as backup
        conversion_successful = False

        if mol_block:
            # Use MOL block
            mol_file = temp_input + ".mol"
            try:
                with open(mol_file, 'w') as f:
                    f.write(mol_block)

                conversion_successful = self._run_openbabel_conversion(
                    mol_file, temp_pdbqt, "mol"
                )
            except Exception as e:
                print(f"Error using MOL block for {mol_id}: {e}")

        # If MOL conversion failed and SMILES is available, try SMILES
        if not conversion_successful and smiles:
            smi_file = temp_input + ".smi"
            try:
                with open(smi_file, 'w') as f:
                    f.write(f"{smiles} {mol_id}\n")

                conversion_successful = self._run_openbabel_conversion(
                    smi_file, temp_pdbqt, "smi", add_3d=True
                )
            except Exception as e:
                print(f"Error using SMILES for {mol_id}: {e}")

        # Copy result to final location if successful
        if conversion_successful and os.path.exists(temp_pdbqt) and os.path.getsize(temp_pdbqt) > 0:
            try:
                with open(temp_pdbqt, 'r') as src, open(output_path, 'w') as dst:
                    dst.write(src.read())
                return True
            except Exception as e:
                print(f"Error copying result for {mol_id}: {e}")

        return False

    def _run_openbabel_conversion(self, input_file, output_file, input_format, add_3d=False):
        """Run OpenBabel conversion with timeout"""
        cmd = ["obabel", input_file, "-O", output_file]

        if add_3d:
            cmd.extend(["--gen3d"])

        # Add hydrogen atoms (important for docking)
        cmd.extend(["-h"])

        try:
            result = subprocess.run(
                cmd, 
                timeout=self.timeout, 
                capture_output=True, 
                text=True
            )

            if result.returncode != 0:
                return False

            return os.path.exists(output_file) and os.path.getsize(output_file) > 0

        except subprocess.TimeoutExpired:
            return False
        except Exception:
            return False

    def convert_library(self, input_file, output_dir=None):
        """Convert molecular library to PDBQT format"""
        # Detect file format and read molecules
        file_format = self.detect_file_format(input_file)

        print(f"Detecting file format: {file_format}")
        print(f"Reading molecules from: {input_file}")

        if file_format == 'smiles':
            molecules = self.read_smiles_file(input_file)
        elif file_format == 'csv':
            molecules = self.read_csv_tsv_file(input_file, delimiter=',')
        elif file_format == 'tsv':
            molecules = self.read_csv_tsv_file(input_file, delimiter='\t')
        elif file_format == 'sdf':
            molecules = self.read_sdf_file(input_file)
        elif file_format == 'mol':
            molecules = self.read_mol_file(input_file)
        else:
            raise ValueError(f"Unsupported format: {file_format}")

        if not molecules:
            print("No molecules found in input file!")
            return

        self.stats['total'] = len(molecules)
        print(f"Found {self.stats['total']} molecules")

        # Setup output directory
        if len(molecules) > 1:
            if output_dir is None:
                output_dir = "ligand"
            os.makedirs(output_dir, exist_ok=True)
            print(f"Output directory: {output_dir}")

        # Convert molecules
        start_time = time.time()

        with tempfile.TemporaryDirectory() as tmpdir:
            for i, molecule in enumerate(molecules, 1):
                # Generate output filename
                if len(molecules) == 1:
                    if output_dir:
                        output_path = os.path.join(output_dir, "lig_0001.pdbqt")
                    else:
                        output_path = f"{Path(input_file).stem}.pdbqt"
                else:
                    output_path = os.path.join(output_dir, f"lig_{i:04d}.pdbqt")

                # Convert molecule
                success = self.convert_molecule_to_pdbqt(molecule, output_path, tmpdir)

                if success:
                    self.stats['converted'] += 1
                else:
                    self.stats['skipped'] += 1

                # Progress update
                processed = self.stats['converted'] + self.stats['skipped']
                elapsed = time.time() - start_time
                rate = processed / elapsed if elapsed > 0 else 0
                eta = (self.stats['total'] - processed) / rate if rate > 0 else 0

                print(f"\rProcessed: {processed}/{self.stats['total']} "
                      f"(converted: {self.stats['converted']}, "
                      f"skipped: {self.stats['skipped']}) "
                      f"Rate: {rate:.1f} mol/s ETA: {eta:.0f}s", end='')

        print("\n" + "="*50)
        print("CONVERSION SUMMARY")
        print("="*50)
        print(f"Total molecules processed: {self.stats['total']}")
        print(f"Successfully converted: {self.stats['converted']}")
        print(f"Skipped (errors/timeouts): {self.stats['skipped']}")
        print(f"Success rate: {(self.stats['converted']/self.stats['total']*100):.1f}%")
        print(f"Total time: {time.time() - start_time:.1f} seconds")

        if self.stats['converted'] > 0:
            if len(molecules) > 1:
                print(f"Output files saved in: {output_dir}/")
            else:
                print(f"Output file: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert molecular libraries to PDBQT format for AutoDock Vina",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python pdbqt_converter.py input.smi
  python pdbqt_converter.py molecules.csv -o output_folder -t 30
  python pdbqt_converter.py library.sdf --timeout 15
        """
    )

    parser.add_argument("input", help="Input molecular file (.smi, .csv, .tsv, .sdf, .mol)")
    parser.add_argument("-o", "--output", help="Output directory (default: 'ligand' for multiple molecules)")
    parser.add_argument("-t", "--timeout", type=int, default=TIMEOUT_SECS,
                       help=f"Timeout for each conversion in seconds (default: {TIMEOUT_SECS})")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found!")
        sys.exit(1)

    try:
        converter = PDBQTConverter(timeout=args.timeout)
        converter.convert_library(args.input, args.output)
    except Exception as e:
        print(f"Error during conversion: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
