import subprocess
import sys
import os
import tempfile

TIMEOUT_SECS = 20  # adjust as preferred

def convert_smiles(smi, mid, mol2file, tmpdir, timeout=TIMEOUT_SECS):
    smi_path = os.path.join(tmpdir, "single.smi")
    mol2_path = os.path.join(tmpdir, "single.mol2")
    with open(smi_path, "w") as f:
        f.write(f"{smi} {mid}\n")
    try:
        result = subprocess.run([
            "obabel", smi_path, "-O", mol2_path, "--gen3d", "-h"
        ], timeout=timeout, capture_output=True, text=True)
        if result.returncode != 0:
            return False  # Open Babel error
        if not os.path.exists(mol2_path) or os.path.getsize(mol2_path) == 0:
            return False
        # Append the single-molecule MOL2 to the output file
        with open(mol2_path) as mol2_in, open(mol2file, "a") as mol2_out:
            mol2_out.write(mol2_in.read())
        return True
    except subprocess.TimeoutExpired:
        # Open Babel took too long—assume stall and skip
        return False

def main(smi_file, mol2_out):
    import shutil
    if os.path.exists(mol2_out):
        os.remove(mol2_out)
    with open(smi_file) as f:
        lines = [l.strip() for l in f if l.strip()]
    total = len(lines)
    count = 0
    skipped = 0
    with tempfile.TemporaryDirectory() as tmpdir:
        for line in lines:
            if " " in line:
                smi, mid = line.split(None, 1)
            else:
                smi, mid = line, "NOID"
            ok = convert_smiles(smi, mid, mol2_out, tmpdir)
            if ok:
                count += 1
            else:
                skipped += 1
            print(f"Processed: {count+skipped}/{total} (kept: {count}, skipped: {skipped})", end='\r')
    print(f"\nDone. Total processed: {count} molecules written, {skipped} skipped.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python babel_smart_convert.py input.smi output.mol2")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

