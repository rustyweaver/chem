#!/bin/env python3

# Extracts some/all of the geometries from an ensemble xyz file
# into another xyz file.
#
# Usage:
#   extracet_geometries.py <input.xyz> <output.xyz> index index index ...
#
# The input file contains a set of .xyz geometries. The indices are
# zero-based.

import sys

def extract_xyz(input_path, output_path, indices):
  with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
    frame_idx = 0
    while True:
      line = fin.readline()
      if not line:
        break
      n_atoms = int(line.strip())
      frame_lines = [line]
      for _ in range(n_atoms + 1):
        frame_lines.append(fin.readline())
      if frame_idx in indices:
        fout.writelines(frame_lines)
      frame_idx += 1

if __name__ == '__main__':
  infile = sys.argv[1]
  outfile = sys.argv[2]
  target_indices = set(int(x) for x in sys.argv[3:])
  extract_xyz(infile, outfile, target_indices)
