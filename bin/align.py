#!/usr/bin/env python3

# Align molecule to be in-plane and centered at the origin.
# This is nasty because you need to know the atom indices of the
# molecules that form the "center", typically a ring. To find the
# atom indices, you can use display.py in the same directory, which
# will create a .html file that displays atom indices in the browser.

import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
import numpy as np

# import numpy as np
# from rdkit import Chem
# from rdkit.Chem import rdDetermineBonds
# import py3Dmol

# 1. Load the molecule from your .xyz file
xyz_filename = "your_trianguleneum_dye.xyz" # Replace with your filename
mol = Chem.MolFromXYZFile(xyz_filename)


# These are the indices of the non-hydrogen atoms that form the
# skeleton of the triangulenium system. But this could easily change if
# the atoms are generated again with smiles or if the .xyz file is
# in any way manipulated.
# These are the two atoms that will sit on the Y axis. The line between
# them is approximately a line of symmetry of the triangulenium skeleton.

MOLECULE_INFO = {
  'nbi-327' : {
    'SKELETON_ATOMS' : [22, 23, 24, 25, 26, 27, 28, 29, 34, 35, 30, 31, 32, 33, 7, 8, 9, 21, 20, 19, 10, 11],
    'ORIGIN_ATOM' : 28,
    'Y_AXIS_ATOM' : 20,
  },
  'nbi-526' : {
    'SKELETON_ATOMS' : [23, 24, 25, 26, 27, 11, 10, 9, 8, 7, 6, 15, 16, 17, 18, 19, 20, 22, 14, 21, 12, 13],
    'ORIGIN_ATOM' : 13,
    'Y_AXIS_ATOM' : 9,
  },
  'htota' : {
    'SKELETON_ATOMS' : [1, 2, 3, 5, 22, 23, 4, 12, 11, 6, 13, 14, 16, 10, 9, 8, 7, 17, 18, 19, 20, 21],
    'ORIGIN_ATOM' : 6,
    'Y_AXIS_ATOM' : 14,
  },
  'tota' : {
    'SKELETON_ATOMS' : [9, 8, 7, 6, 11, 10, 5, 4, 12, 13, 3, 2, 1, 0, 21, 15, 14, 16, 17, 18, 19, 20],
    'ORIGIN_ATOM' : 13,
    'Y_AXIS_ATOM' : 17,
  },
  'daota' : {
    'SKELETON_ATOMS' : [7, 8, 9, 10, 11, 6, 20, 19, 14, 13, 18, 17, 16, 15, 21, 0, 12, 1, 2, 3, 4, 5],
    'ORIGIN_ATOM' : 13,
    'Y_AXIS_ATOM' : 2,
  },
  'tata' : {
    'SKELETON_ATOMS' : [7, 8, 9, 10, 11, 6, 20, 19, 14, 13, 18, 17, 16, 15, 21, 0, 12, 1, 2, 3, 4, 5],
    'ORIGIN_ATOM' : 13,
    'Y_AXIS_ATOM' : 2,
  },
  'adota' : {
    'SKELETON_ATOMS' : [7, 8, 9, 10, 11, 6, 20, 19, 14, 13, 18, 17, 16, 15, 21, 0, 12, 1, 2, 3, 4, 5],
    'ORIGIN_ATOM' : 13,
    'Y_AXIS_ATOM' : 2,
  },
  'aadota' : {
    'SKELETON_ATOMS' : [22, 24, 18, 19, 20, 21, 17, 16, 7, 6, 15, 13, 12, 8, 9, 10, 5, 11, 1, 2, 3, 4],
    'ORIGIN_ATOM' : 6,
    'Y_AXIS_ATOM' : 13,
  },
  'atota' : {
    'SKELETON_ATOMS' : [13, 15, 16, 7, 8, 12, 17, 18, 19, 5, 24, 22, 21, 20, 4, 3, 5, 2, 1, 11, 10, 9],
    'ORIGIN_ATOM' : 6,
    'Y_AXIS_ATOM' : 13,
  },
  'guanine' : {
    'SKELETON_ATOMS' : [13, 15, 16, 7, 8, 12, 17, 18, 19, 5, 24, 22, 21, 20, 4, 3, 5, 2, 1, 11, 10, 9],
    'ORIGIN_ATOM' : 6,
    'Y_AXIS_ATOM' : 13,
  },
}


# Translates the molecule by moving the center of the skeleton
# to the origin.
def translate_to_origin(mol, origin_atom_index):
  mol = Chem.Mol(mol)
  conf = mol.GetConformer()
  coords = conf.GetPositions()
  origin_coord = coords[origin_atom_index]
  new_coords = coords - origin_coord
  conf.SetPositions(new_coords)
  return mol


# Rotates the molecule so that the (approximately) planar skeleton
# sits in the x-y plane.
def align_to_xy(mol, skeleton_indices):
  mol = Chem.Mol(mol)
  conf = mol.GetConformer()
  coords = conf.GetPositions()
  _, _, vh = np.linalg.svd(coords[skeleton_indices])
  normal = vh[2]
  z_axis = np.array([0.0, 0.0, 1.0])
  v = np.cross(normal, z_axis)
  s = np.linalg.norm(v)
  c = np.dot(normal, z_axis)
  vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
  R1 = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s**2)) if s > 1e-6 else np.eye(3)
  conf.SetPositions(coords @ R1.T)
  return mol


# Rotates the molecule about the z axis until the line of symmetry
# (defined by two atom_indices) is coincident with the y axis.
# The origin_atom_index should identify an atom that is already at
# the orign.
def rotate_to_y(mol, origin_atom_index, y_axis_atom_index):
  mol = Chem.Mol(mol)
  conf = mol.GetConformer()
  coords = conf.GetPositions()
  assert np.allclose(coords[origin_atom_index], [0.0, 0.0, 0.0], atol=1e-5)
  diff = coords[y_axis_atom_index] - coords[origin_atom_index]
  phi = np.arctan2(diff[1], diff[0])
  theta = (np.pi / 2) - phi
  cos_a, sin_a = np.cos(theta), np.sin(theta)
  R = np.array([[cos_a, -sin_a, 0], [sin_a, cos_a, 0], [0, 0, 1]])
  conf.SetPositions(coords @ R.T)
  return mol


def fatal(s):
  print("ERROR: %s" % s, file=sys.stderr)
  sys.exit(1)

def main():
  if len(sys.argv) < 4:
    fatal("usage: rotate.py input.xyz system output.xy")

  input_xyz = sys.argv[1]
  system = sys.argv[2]
  output_xyz = sys.argv[3]

  if not system in MOLECULE_INFO:
    fatal("Unknown system %s must be one of: %s" % (system, " ".join(MOLECULE_INFO.keys())))

  if not Path(input_xyz).exists():
    fatal("Input file %s not found" % input_xyz)
  print('Reading', input_xyz)
  input_mol = Chem.MolFromXYZFile(input_xyz)
  if not input_mol:
    fatal("Unable to read XYZ file %s" % input_xyz)
  mol_info = MOLECULE_INFO[system]
  skeleton_atoms = mol_info['SKELETON_ATOMS']
  origin_atom = mol_info['ORIGIN_ATOM']
  y_axis_atom = mol_info['Y_AXIS_ATOM']
  output_mol = rotate_to_y(
    align_to_xy(translate_to_origin(input_mol, origin_atom),
                skeleton_atoms), origin_atom, y_axis_atom)
  print('Writing', output_xyz)
  Chem.MolToXYZFile(output_mol, output_xyz)
  
if __name__ == "__main__":
    main()
