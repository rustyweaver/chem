#!/usr/bin/env python3
import sys
import re

def parse_hly_charges(text):
    """
    Return the raw text of the LAST ESP/HLY charges block,
    exactly as it appears in the Gaussian output.
    """

    # Capture EVERYTHING from "ESP charges:" through "Sum of ESP charges = ..."
    matches = re.findall(
        r"(ESP charges:.*?Sum of ESP charges[^\n]*)",
        text,
        flags=re.S
    )

    if not matches:
        return ""

    # Return the last matching block, unchanged
    return matches[-1]
