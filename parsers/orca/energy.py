#!/usr/bin/env python3
import sys
import re

"""
-------------------------   --------------------
FINAL SINGLE POINT ENERGY       -93.351294561124
-------------------------   --------------------
"""

expr = (
    r"^\s*FINAL SINGLE POINT ENERGY\s+(\S+)\s*$"
)


def parse_energy(text):
    """
    Return final single point energy from orca output.
    """

    matches = re.findall(expr, text, re.M)
    if not matches:
        return None
    else:
        return float(matches[-1])
