#!/usr/bin/env python3
import sys
import re

"""
Sum of electronic and thermal Free Energies=         -115.682144
"""

gibbs_energy_expr = (
    r"^\s*Sum of electronic and thermal Free Energies=\s+(\S+)\s*$"
)

def parse_gibbs_energy(text):
    """
    Return final single point energy from orca output.
    """

    matches = re.findall(gibbs_energy_expr, text, re.M)
    if not matches:
        return None
    else:
        return float(matches[-1])

"""
Sum of electronic and thermal Enthalpies=           -115.655078
"""

enthalpy_expr = (
    r"^\s*Sum of electronic and thermal Enthalpies=\s+(\S+)\s*$"
)

def parse_enthalpy(text):
    """
    Return final single point energy from orca output.
    """

    matches = re.findall(enthalpy_expr, text, re.M)
    if not matches:
        return None
    else:
        return float(matches[-1])
