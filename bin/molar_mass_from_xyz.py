#!/usr/bin/env python3

import sys
from ase.io import read

def main():
  if len(sys.argv) < 2:
    print("Usage: molar_mass_from_xyz.py <molecule.xyz>")
    sys.exit(1)

  filename = sys.argv[1]
  atoms = read(filename)
  molar_mass = sum(atoms.get_masses())
  print(f"Total Molar Mass: {molar_mass:.4f} g/mol\n")

  counts = {}
  for atom in atoms:
    sym = atom.symbol
    if sym in counts:
      counts[sym] += 1
    else:
      counts[sym] = 1
  print(counts)

  # print(f"Atom: {atom.symbol}, Position: {atom.position}")

if __name__ == "__main__":
  main()
