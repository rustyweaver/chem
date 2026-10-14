#!/usr/bin/env python3
import sys
import re

"""
Total Dipole Moment    :     -0.496880959      -0.081840842      -0.148613552
                        -----------------------------------------
Magnitude (a.u.)       :      0.525047234
Magnitude (Debye)      :      1.334564010
"""

# The pattern captures 4 lines exactly as requested
expr = (
    r"^Total Dipole Moment.*\n"       # Line 1: Starts with string, matches rest of line + newline
    r".*\n"                           # Line 2: Matches entire line + newline
    r".*\n"                           # Line 3: Matches entire line + newline
    r"Magnitude \(Debye\)\s*:\s*(\S+).*" # Line 4: Capture group, then match rest of line
)

def parse_dipole(text):
    # re.MULTILINE (re.M) is essential so ^ matches the start of EACH line
    vals = re.findall(expr, text, re.M)
    return [float(v) for v in vals] if vals else None
