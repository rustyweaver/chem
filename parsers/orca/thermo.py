#!/usr/bin/env python3
import sys
import re

"""
-------------------
GIBBS FREE ENERGY
-------------------

Total enthalpy                    ...  -1068.86967731 Eh
Total entropy correction          ...     -0.06897906 Eh    -43.28 kcal/mol
-----------------------------------------------------------------------
Final Gibbs free energy         ...  -1068.93865637 Eh
"""

expr_H = r"^\s*Total enthalpy\s+\.\.\.\s+(\S+)\s+Eh\s*$"
expr_TS = r"^\s*Total entropy correction\s+\.\.\.\s+(\S+)\s+Eh"
expr_G = r"^\s*Final Gibbs free energy\s+\.\.\.\s+(\S+)\s+Eh\s*$"


def parse_thermo(text):
    """
    Return (H, -T*S, G) from ORCA output in Eh.
    """
    H = re.findall(expr_H, text, re.M)
    TS = re.findall(expr_TS, text, re.M)
    G = re.findall(expr_G, text, re.M)

    if not (H and TS and G):
        return None

    return {
        "enthalpy": float(H[-1]),
        "entropy_correction": float(TS[-1]),
        "gibbs_energy": float(G[-1]),
    }
