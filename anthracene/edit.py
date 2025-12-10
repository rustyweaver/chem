#!/usr/bin/env python3
#
# Add aromatic bonds to fused rings preserving the atom indices.
# will.rusch@gmail.com

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Chem
import numpy as np

# The smiles string for bis(anthracenyl)benzene.
smiles = "c1(c4c(c5c6ccccc6cc7ccccc75)cccc4)c2ccccc2cc3ccccc31"

# From the smiles derive the smarts pattern
smarts = "[c:100]1(c4c([c:200]5c6ccccc6[c:300]c7ccccc75)cccc4)c2ccccc2[c:400]c3ccccc31"

# Notes:
# If you want to adapt this code to use a molecule that contains
# bis(anthracenuyl)benzene plus some other subsituents, functional
# groups, etc, you need to change 'c' in the pattern below with
# 'c(your functional group)' and you probably also read the smiles
# syntax documents, esp. regarding lower case c (aromatic) and upper
# case C (explicit). We recommend just eliminating all hydrogens here
# because they will be automatically added later.
#
# Create the smarts string. This is similar to smiles but it serves
# a different purpose: pattern matching. It is just like smiles but adds some
# fancy syntax for matching some other molecules with the same core pattern.
#
# Attention: a 'c' in a smarts string matches carbon or carbon with H (C-H),
# or carbon with any R group (C-R).
#
# Attention: we will *label* certain atoms in the pattern so we can track
# then down later and add bonds. We will add bonds at four positions, but
# to visualize that you really need to look at this article / image.
# The labeling syntax is [c:XXX] where XXX is an arbitrary integer label that you
# invent and later use to locate that atom in the specific molecule used.
#
# To the original molecule We add two bonds at positions 1 <-> 2 and 3 <-> 4.
# We call this new create the 'excited conformer'.
#
# For explanation only, we mark the four bonding atoms in the smiles string
# with a #.
#
#   "c#1(c4c(c#5c6#ccccc6c#c7ccccc75)cccc4)c2c#ccc2cc3ccccc31"
#
# In the smarts string, we actually need to use the syntax [c:101] to label
# that specific carbon as 101 (for example). We thus substitute 'c' with '[c:101]'
# again, for example.


# Creates a mapping from the matching label (eg. 100) to the
# corresponding atom in the smiles string.

def get_match_map(pattern, match):
  mapped = {}
  for atom in pattern.GetAtoms():
    amap = atom.GetAtomMapNum()
    if amap > 0:
      # amap equals the custom index from the smarts
      # pattern, eg. [X:1] then amap = 1.
      mapped[amap] = match[atom.GetIdx()]
  return mapped


# Adds a bond between two carbons. It will also unbond an attached
# hydrogen (if present) on each side of the bond.  In order to
# preserve the atom indices, the hydrogens are first unbonded, and
# then removed to a far away location, but remain in molecule.

def connect_atoms_preserve_indices(molH, a, b, bond_type=Chem.rdchem.BondType.SINGLE):
    rw = Chem.RWMol(molH)

    # copy conformer
    conf = molH.GetConformer()
    rw.RemoveAllConformers()
    rw.AddConformer(conf, assignId=True)
    conf = rw.GetConformer()

    def find_H(idx):
        for nbr in rw.GetAtomWithIdx(idx).GetNeighbors():
            if nbr.GetSymbol() == "H":
                return nbr.GetIdx()
        return None

    h1 = find_H(a)
    h2 = find_H(b)

    if h1 is not None:
        rw.RemoveBond(a, h1)
    if h2 is not None:
        rw.RemoveBond(b, h2)

    rw.AddBond(a, b, bond_type)

    # move hydrogens away
    far = np.array([999.0,999.0,999.0])
    if h1 is not None:
        conf.SetAtomPosition(h1, far)
    if h2 is not None:
        conf.SetAtomPosition(h2, far + np.array([1.0,0,0]))

    mol_new = rw.GetMol()

    # do not sanitize fully (avoids kekulization)
    Chem.SanitizeMol(mol_new, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ADJUSTHS)

    return mol_new


def write_xyz(mol, filename):
    xyz = Chem.MolToXYZBlock(mol)
    with open(filename, "w") as f:
        f.write(xyz)


def write_mol(mol, filename):
    molblock = Chem.MolToMolBlock(mol, kekulize=False)
    with open(filename, "w") as f:
        f.write(molblock)


mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

#  generate 3d coords before connecting atoms
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol)

pattern = Chem.MolFromSmarts(smarts)
match = mol.GetSubstructMatch(pattern)
assert match, ( 'The smiles string "%s" does not match the ' \
                'bis(anthracenyl)benzene pattern' % smiles )

match_map = get_match_map(pattern, match)
newmol = connect_atoms_preserve_indices(mol, match_map[100], match_map[200])
newmol = connect_atoms_preserve_indices(newmol, match_map[300], match_map[400])

write_xyz(mol, "mol.xyz")
write_mol(mol, "mol.mol")

write_xyz(newmol, "newmol.xyz")
write_mol(newmol, "newmol.mol")
