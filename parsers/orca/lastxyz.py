#!/usr/bin/env python3

import argparse
import json
from filename import parse_filename
from absorption import parse_absorption_spectrum
from energy import parse_energy
from thermo import parse_thermo
from dipole import parse_dipole
from bond import parse_bond
from geometry import get_final_geometry

def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("files", nargs="+")
  tasks = set()
  args = ap.parse_args()

  records = []
  for f in args.files:
    with open(f) as h:
      text = h.read()
      print(get_final_geometry(text))

if __name__ == "__main__":
  main()
