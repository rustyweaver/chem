from rdkit import Chem
from rdkit.Chem import AllChem, Draw

# bis(anthracenyl)benzene
smiles = "c1(c4c(c5c6ccccc6cc7ccccc75)cccc4)c2ccccc2cc3ccccc31"

# 2D drawing (NO hydrogens)
mol2d = Chem.MolFromSmiles(smiles)
Draw.MolToFile(mol2d, "molecule.png", size=(600,600))

# 3D version (WITH hydrogens)
mol3d = Chem.MolFromSmiles(smiles)
mol3d = Chem.AddHs(mol3d)
AllChem.EmbedMolecule(mol3d, AllChem.ETKDG())
Chem.MolToMolFile(mol3d, "molecule.mol")
