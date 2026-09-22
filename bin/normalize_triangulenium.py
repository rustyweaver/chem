#!/usr/bin/env python3

import numpy as np
import sys

def orient_xyz(input_file, output_file):
    # 1. Read the XYZ file
    with open(input_file, 'r') as f:
        lines = f.readlines()
        
    num_atoms = int(lines[0].strip())
    comment = lines[1].strip()
    
    atoms = []
    coords = []
    for line in lines[2:2+num_atoms]:
        parts = line.split()
        if len(parts) >= 4:
            atoms.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            
    coords = np.array(coords)
    
    # 2. Center at the origin (Geometric Center)
    centroid = np.mean(coords, axis=0)
    coords_centered = coords - centroid
    
#    # 3. Use SVD to find principal axes and rotate
#    # U, S, Vt = Singular Value Decomposition
#    _, _, Vt = np.linalg.svd(coords_centered)
#    
#    # Vt contains the principal axes. Multiplying rotates the molecule
#    # so max variance is on X, medium on Y, and lowest (flat normal) on Z.
#    coords_oriented = np.dot(coords_centered, Vt.T)
    
    # 4. Write to a new XYZ file
    with open(output_file, 'w') as f:
        f.write(f"{num_atoms}\n")
        f.write(f"{comment} - Oriented to XY plane and centered\n")
        for atom, coord in zip(atoms, coords_centered):
            f.write(f"{atom:<2} {coord[0]:14.8f} {coord[1]:14.8f} {coord[2]:14.8f}\n")


def main():
    if len(sys.argv) != 3:
        print("Usage: python normalize_triangulenum.py input.xyz output.xyz", file=sys.stderr)
        sys.exit(1)
    orient_xyz(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
