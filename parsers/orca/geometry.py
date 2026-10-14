#!/usr/bin/env python3

import re

# Block 1: The convergence marker
CONVERGENCE_HEADER = """\
                    ***********************HURRAY********************
                    *** THE OPTIMIZATION HAS CONVERGED     ***
                    *************************************************
"""

# Block 2: The coordinate section marker
COORDS_HEADER = """
*******************************************************
*** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***
*** (AFTER  .* CYCLES)               ***
*******************************************************
"""

def extract_exact_block(raw_block):
    """Converts a raw pasted block into a robust regex pattern while allowing .* wildcards."""
    lines = [l for l in raw_block.splitlines() if l.strip()]
    processed = []
    for line in lines:
        # 1. Escape the whole line first to protect (, ), *, etc.
        escaped = re.escape(line.strip())
        # 2. Re-enable the wildcard by un-escaping any \.\* back to .*
        flexible_wildcard = escaped.replace(r'\.\*', '.*')
        # 3. Turn escaped spaces into flexible whitespace
        flexible_space = re.sub(r'(\\ )+', r'\\s+', flexible_wildcard)
        processed.append(flexible_space)
    return r"\s+".join(processed)

def get_final_geometry(s):
    # 1. Find everything after the 'HURRAY' convergence block
    conv_pattern = extract_exact_block(CONVERGENCE_HEADER) + r"\s+(.*)"
    conv_match = re.search(conv_pattern, s, re.DOTALL)

    if not conv_match:
        raise ValueError("Convergence 'HURRAY' block not found.")

    s1 = conv_match.group(1)

    # 2. Match the header and capture the coordinate block
    coords_pattern = extract_exact_block(COORDS_HEADER) + r"\s+(.*?)(?:\n\s*\n|$)"
    coords_match = re.search(coords_pattern, s1, re.DOTALL)

    if not coords_match:
        raise ValueError("Stationary Point block not found.")

    # 3. FIX: Split lines and skip the 3 header lines
    raw_lines = coords_match.group(1).strip().splitlines()
    clean_coords = "\n".join(raw_lines[3:]) 
    
    return clean_coords.strip()
