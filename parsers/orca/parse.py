#!/usr/bin/env python3

import argparse
import csv
import json
import sys
from filename import parse_filename
from absorption import parse_absorption_spectrum
from energy import parse_energy
from thermo import parse_thermo
from dipole import parse_dipole
from bond import parse_scanned_bond
from bond import parse_bond_length
from bond import parse_bond_angle
from geometry import get_final_geometry

def parse_bond_indices(val):
  try:
    i, j = map(int, val.split(","))
    if i > 0 and j > 0:
      return (i, j)
  except ValueError:
    pass
  raise argparse.ArgumentTypeError("Bond length requires two positive integers separated by a comma (e.g. 1,2)")

def parse_angle_indices(val):
  print('getting angle indices')
  try:
    i, j, k = map(int, val.split(","))
    if i > 0 and j > 0 and k > 0:
      return (i, j, k)
  except ValueError:
    pass
  raise argparse.ArgumentTypeError("Bond angle requires three positive integers separated by a comma (e.g. 1,2,3)")

def flatten_record(record):
  list_keys = [k for k, v in record.items() if isinstance(v, list) and v and isinstance(v[0], dict)]
  if list_keys:
    list_key = list_keys[0]
    base_meta = {k: v for k, v in record.items() if k != list_key}
    return [{**base_meta, **item} for item in record[list_key]]
  return [record]

def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("files", nargs="+")
  ap.add_argument("--format", choices=["json", "csv"], default="csv",
                  help="output format (json or csv)")
  ap.add_argument("--columns", type=lambda s: [c.strip() for c in s.split(",")],
                  help="comma-separated list of columns to include")
  ap.add_argument("--absorption", action="store_true",
                  help="whether to capture the absorption spectrum")
  ap.add_argument("--energy", action="store_true",
                  help="final single point energy")
  ap.add_argument("--thermo", action="store_true",
                  help="final single point thermo")
  ap.add_argument("--dipole", action="store_true",
                  help="dipole moment")
  ap.add_argument("--scanned_bond", action="store_true",
                  help="bond length when scanning")
  ap.add_argument("--geometry", action="store_true",
                  help="final converged geometry in Angstrom")
  ap.add_argument("--bond_length", type=parse_bond_indices, metavar="i,j",
                  help="distance between two atoms i,j")
  ap.add_argument("--bond_angle", type=parse_angle_indices, metavar="i,j,k",
                  help="bond angle between three atoms i,j,k")

  tasks = {}
  args = ap.parse_args()

  if args.absorption:
    tasks["absorption"] = lambda text: parse_absorption_spectrum(text)
  if args.energy:
    tasks["energy"] = lambda text: parse_energy(text)
  if args.thermo:
    tasks["thermo"] = lambda text: parse_thermo(text)
  if args.dipole:
    tasks["dipole"] = lambda text: parse_dipole(text)
  if args.scanned_bond:
    tasks["scanned_bond"] = lambda text: parse_scanned_bond(text)
  if args.geometry:
    tasks["geometry"] = lambda text: get_final_geometry(text)
  if args.bond_length:
    tasks["bond_length"] = lambda text: parse_bond_length(text, *args.bond_length)
  if args.bond_angle:
    tasks["bond_angle"] = lambda text: parse_bond_angle(text, *args.bond_angle)

  writer = None
  for f in args.files:
    record = parse_filename(f)
    with open(f) as h:
      text = h.read()
    for task_name, parser in tasks.items():
      record[task_name] = parser(text)

    flat_rows = flatten_record(record)
    if args.columns:
      flat_rows = [{k: row[k] for k in args.columns if k in row} for row in flat_rows]

    if args.format == "json":
      for row in flat_rows:
        print(json.dumps(row))
    else:
      if writer is None and flat_rows:
        fieldnames = args.columns if args.columns else list(flat_rows[0].keys())
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
        writer.writeheader()
      if writer and flat_rows:
        writer.writerows(flat_rows)

if __name__ == "__main__":
  main()
