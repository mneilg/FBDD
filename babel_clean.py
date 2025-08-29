import subprocess
import os
import tempfile
import sys

TIMEOUT_SECS = 20 # adjust as preferred

def process_sdf_molecules(sdf_file, sdf_out, timeout=TIMEOUT_SECS):
    if os.path.exists(sdf_out):
        os.remove(sdf_out)

    with open(sdf_file) as f:
        content = f.read()
        molecules = content.split('$$$$\n')   # SDF delimiter

    total = len([m for m in molecules if m.strip()])
    count = 0
    skipped = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        for idx, molecule in enumerate(molecules):
            if not molecule.strip():
                continue
            molecule_block = molecule if molecule.endswith('$$$$\n') else molecule + '\n$$$$\n'
            molecule_path = os.path.join(tmpdir, f"mol_{idx}.sdf")
            converted_path = os.path.join(tmpdir, f"conv_{idx}.sdf")
            with open(molecule_path, "w") as mf:
                mf.write(molecule_block)
            try:
                result = subprocess.run(
                    [
                        "obabel", molecule_path, "-O", converted_path, "--gen3d"
                    ],
                    timeout=timeout, capture_output=True, text=True
                )
                if result.returncode != 0 or not os.path.exists(converted_path) or os.path.getsize(converted_path) == 0:
                    skipped += 1
                    continue
                with open(converted_path) as conv_mf, open(sdf_out, "a") as final_out:
                    conv_content = conv_mf.read()
                    final_out.write(conv_content)
                count += 1
            except subprocess.TimeoutExpired:
                skipped += 1
                continue
            print(f"Processed: {count+skipped}/{total} (kept: {count}, skipped: {skipped})", end='\r')
    print(f"\nDone. Total processed: {count} molecules written, {skipped} skipped.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python sdf_clean_convert.py input.sdf output.sdf")
        sys.exit(1)
    process_sdf_molecules(sys.argv[1], sys.argv[2])


