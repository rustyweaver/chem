#!/usr/bin/env python3

import argparse
import sys
import py3Dmol
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

parser = argparse.ArgumentParser()
parser.add_argument("inputfile")
parser.add_argument("outputfile")
parser.add_argument("--charge", type=int, default=0)
parser.add_argument("--show_atom_labels", action="store_true")
parser.add_argument("--show_axes", action="store_true")
parser.add_argument(
  "--tdm",
  type=lambda s: [float(x) for x in s.split(",")],
  action="append",
  default=[],
)

args = parser.parse_args()
for entry in args.tdm:
  assert len(entry) == 4, "Each tdm argument must be four comma-separated floats"

mol = Chem.MolFromXYZFile(args.inputfile)
rdDetermineBonds.DetermineBonds(mol, charge=args.charge)
xyz_block = Chem.MolToXYZBlock(mol)
view = py3Dmol.view(width=800, height=800)
view.addModel(xyz_block, "xyz")

# Label the atoms with the index of the atom in the .xyz file.
# Indices are zero-based.
if args.show_atom_labels:
  conf = mol.GetConformer()
  for atom in mol.GetAtoms():
    idx = atom.GetIdx()
    pos = conf.GetAtomPosition(idx)
    view.addLabel(
      str(idx),
      {
        "position": {"x": pos.x, "y": pos.y, "z": pos.z},
        "showBackground": False,
        "fontSize": 10,
      },
    )


view.setStyle({'not': {'elem': 'H'}}, {'stick': {'radius': 0.15}, 'sphere': {'scale': 0.3}})


if args.show_axes:

  axis_length = 10.0
  radius = 0.04
  
  # X-axis (red arrow)
  view.addArrow({
      "start": {"x": 0, "y": 0, "z": 0},
      "end": {"x": axis_length, "y": 0, "z": 0},
      "radius": 0.04,
      "color": "red",
      # distance from the end of the arrow where where the pointy
      # part begins.
      "midpos": -0.4
   })
  view.addLabel(
      "x",
      {
          "position": {"x": axis_length + 0.4, "y": 0, "z": 0},
          "showBackground": False,
          "fontColor": "red",
          "fontSize": 14,
      },
   )
  
  # Y-axis (green arrow)
  view.addArrow({
      "start": {"x": 0, "y": 0, "z": 0},
      "end": {"x": 0, "y": axis_length, "z": 0},
      "radius": 0.04,
      "color": "green",
      # distance from the end of the arrow where where the pointy
      # part begins.
      "midpos": -0.4
   })
  view.addLabel(
      "y",
      {
          "position": {"x": 0, "y": axis_length + 0.4, "z": 0},
          "showBackground": False,
          "fontColor": "green",
          "fontSize": 14,
      },
   )
  
  # Z-axis (blue arrow)
  view.addArrow({
      "start": {"x": 0, "y": 0, "z": 0},
      "end": {"x": 0, "y": 0, "z": axis_length},
      "radius": 0.04,
      "color": "blue",
      # distance from the end of the arrow where where the pointy
      # part begins.
      "midpos": -0.4
   })
  view.addLabel(
      "z",
      {
          "position": {"x": 0, "y": 0, "z": axis_length + 0.4},
          "showBackground": False,
          "fontColor": "blue",
          "fontSize": 14,
      },
   )

def mag(x, y, z):
  return (x**2 + y**2 + z**2) ** 0.5

# Display transition dipole moments if given
colors = ["purple", "orange", "cyan", "magenta", "yellow"]
for i, entry in enumerate(args.tdm):
  osc, tx, ty, tz = entry
  tmag = mag(tx, ty, tz)
  if tmag > 0:
    tx /= tmag
    ty /= tmag
    tz /= tmag
  tx *= osc * 30
  ty *= osc * 30
  tz *= osc * 30
  view.addArrow({
    "start": {"x": 0, "y": 0, "z": 0},
    "end": {"x": tx, "y": ty, "z": tz},
    "radius": 0.1,
    "color": colors[i % len(colors)],
    "midpos": -0.4,
  })

view.zoomTo()

print("Writing:", args.outputfile)
with open(args.outputfile, "w") as f:
  f.write(view._make_html())
