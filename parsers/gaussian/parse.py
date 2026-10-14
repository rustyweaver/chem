#!/usr/bin/env python3

import argparse
import json
from filename import parse_filename
from mulliken import parse_mulliken_charges
from esp import parse_esp_charges
from dipole import parse_dipole_moment
from flatten import flatten_all_by_charge
from energy import parse_gibbs_energy
from energy import parse_enthalpy

def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("files", nargs="+")
  ap.add_argument("--charges", action="append", default=[],
                  help="Which charge models to extract: mulliken, esp, all")
  ap.add_argument("--dipole", action="store_true", help="Extract the dipole moment")
  ap.add_argument("--gibbs_energy", action="store_true", help="Extract the Gibbs energy")
  ap.add_argument("--enthalpy", action="store_true", help="Extract the Enthalpy")

  args = ap.parse_args()
  dispatch = {
    "mulliken": parse_mulliken_charges,
    "esp": parse_esp_charges,
    "dipole_moment": parse_dipole_moment,
    "gibbs_energy": parse_gibbs_energy,
    "enthalpy": parse_enthalpy
  }

  requested = set(args.charges)
  if "all" in requested:
    requested.update({"mulliken", "esp"})
  if args.dipole: requested.add("dipole_moment")
  if args.gibbs_energy: requested.add("gibbs_energy")
  if args.enthalpy: requested.add("enthalpy")
  
  for f in args.files:
    record = parse_filename(f)
    with open(f) as h:
      text = h.read()
    for key in (requested & dispatch.keys()):
      record[key] = dispatch[key](text)
    print(json.dumps(record))

if __name__ == "__main__":
  main()
