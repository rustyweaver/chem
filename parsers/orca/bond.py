#!/usr/bin/env python3

import math
import re
from geometry import get_final_geometry

"""
         *               RELAXED SURFACE SCAN STEP   1               *
         *                                                           *
         *                 Bond (  1,   0)  :   0.50000000           *
         *************************************************************
"""

expr = (
    r"^\s*\*\s*RELAXED SURFACE SCAN STEP\s+\d+\.?\d*.*\n"
    r".*\n"
    r"\s*\*\s*Bond\s*\(.*\)\s*:\s*(\d+\.?\d*).*"
)

def parse_scanned_bond(text):
    vals = re.findall(expr, text, re.M)
    return [float(v) for v in vals] if vals else None


def get_bond_length_from_xyz_formatted_text(text, i, j):
  lines = [line.split() for line in text.strip().splitlines() if line.split()]
  coords = [[float(x) for x in line[1:4]] for line in lines]
  p1, p2 = coords[i - 1], coords[j - 1]
  return math.dist(p1, p2)


def parse_bond_length(text, i, j):
  geom = get_final_geometry(text)
  return get_bond_length_from_xyz_formatted_text(geom, i, j)


def get_bond_angle_from_xyz_formatted_text(text, i, j, k):
  lines = [line.split() for line in text.strip().splitlines() if line.split()]
  coords = [[float(x) for x in line[1:4]] for line in lines]
  v1 = [coords[i - 1][m] - coords[j - 1][m] for m in range(3)]
  v2 = [coords[k - 1][m] - coords[j - 1][m] for m in range(3)]
  dot = sum(a * b for a, b in zip(v1, v2))
  mag1 = math.sqrt(sum(a**2 for a in v1))
  mag2 = math.sqrt(sum(a**2 for a in v2))
  return math.degrees(math.acos(max(-1.0, min(1.0, dot / (mag1 * mag2)))))


def parse_bond_angle(text, i, j, k):
  '''Returns the bond angle of atoms i, j, k in degrees.'''
  geom = get_final_geometry(text)
  return get_bond_angle_from_xyz_formatted_text(geom, i, j, k)
