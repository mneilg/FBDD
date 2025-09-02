import subprocess
import os
import tempfile
import sys
import re

TIMEOUT_SECS = 20  # adjust as preferred

def extract_sdf_fields(sdf_block):
    # Find all fields of the form >  <Tag>
    fields = {}
    # Regular expression for SDF field blocks
    field_blocks = re.findall(r'> *<([^>]+)>[\r\n]+(.*?)[\r\n]{2,}', sdf_block, re.DOTALL)
    for tag, val in field_blocks:
        fields[tag.strip()] = val.strip()
    return fields

def append_sdf_fields(sdf_block, fields):
    # Remove old SDF fields first, then append fresh ones at the end (before $$$$)
    sdf_block = re.sub(r'> *<[^>]+>[\r\n]+.*?[\r\n]{2,}', '', sdf_block, flags=re.DOTALL)
    sdf_block = sdf_block.rstrip()
    for tag, val in fields.items():
        sdf_block += f"\n>  <{tag}>\n{val}\n"
    sdf_block += "\n$$$$\n"
    return sdf_block

def process_sdf_molecules(sdf_file, sdf_out, timeout=TIMEOUT_SECS):
    if os.path.exists(sdf_out):
        os.remove(sdf_out)

    with open(sdf_file) as f:
        content = f.read()

    molecules = [m for m in content.split('$$$$') if m.strip()]
    total = len(molecules)
    count = 0
    skipped = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        for idx, molecule in enumerate(molecules):
            molecule_block = molecule.strip()
            if not molecule_block: continue
            # Ensure original SDF formatting
            if not molecule_block.endswith('\n'):
                molecule_block += '\n'
            molecule_block += '$$$$\n'
            # Extract fields before conversion
            sdf_fields = extract_sdf_fields(molecule_block)

            molecule_path = os.path.join(tmpdir, f"mol_{idx}.sdf")
            converted_path = os.path.join(tmpdir, f"conv_{idx}.sdf")
            with open(molecule_path, "w") as mf:
                mf.write(molecule_block)

            try:
                result = subprocess.run([
                        "obabel", molecule_path, "-O", converted_path, "--gen3d"
                    ],
                    timeout=timeout, capture_output=True, text=True
                )
                if result.returncode != 0 or not os.path.exists(converted_path) or os.path.getsize(converted_path) == 0:
                    skipped += 1
                    continue
                # Read converted block (typically Open Babel writes only the structure, may lose fields)
                with open(converted_path) as conv_mf:
                    conv_content = conv_mf.read()
                # Add back property fields to output molecule
                final_block = append_sdf_fields(conv_content.split('$$$$')[0], sdf_fields)
                with open(sdf_out, "a") as final_out:
                    final_out.write(final_block)
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

