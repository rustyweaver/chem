#!/usr/bin/env python3
import sys
import re

def parse_outfile(outfile, cutoff, max_conformers):
    """Return list of conformer indices up to cutoff% cumulative population."""
    conformers = []
    pattern = re.compile(r'^\s*(\d+)\s+([\d\.\-]+)\s+\d+\s+([\d\.]+)\s+([\d\.]+)')
    with open(outfile) as f:
        for line in f:
            m = pattern.match(line)
            if m:
                idx = int(m.group(1))
                cumul = float(m.group(4))
                conformers.append((idx, cumul))

    selected = []
    for idx, cumul in conformers:
      selected.append(idx)
      if cumul >= cutoff or len(selected) > max_conformers:
        break;

    print(f"{outfile}: selected {len(selected)} conformers up to {cutoff}% cumulative population.")
    return selected

def split_xyz(filename, prefix, selected):
    """Extract and write only the selected conformers from a multi-conformer XYZ file."""
    with open(filename) as f:
        lines = f.readlines()

    i = 0
    conf = 0
    written = set()
    while i < len(lines):
        n_atoms = int(lines[i].strip())
        block = lines[i:i + n_atoms + 2]
        if conf in selected:
            out_name = f"{prefix}-conf-{conf:02d}.xyz"
            with open(out_name, "w") as out:
                out.writelines(block)
            written.add(conf)
        i += n_atoms + 2
        conf += 1

    missing = [c for c in selected if c not in written]
    assert not missing, f"Missing conformers from {filename}: {missing}"

def main():
    if len(sys.argv) != 6:
        print("Usage: python split-goat-xyz.py <input.xyz> <prefix> <out file> <cutoff%> <max-conformers>")
        sys.exit(1)
    xyz_file, prefix, out_file, cutoff, max_conformers = sys.argv[1:6]
    cutoff = float(cutoff)
    max_conformers = int(max_conformers)
    selected = parse_outfile(out_file, cutoff, max_conformers)
    split_xyz(xyz_file, prefix, selected)

if __name__ == "__main__":
    main()
