#!/usr/bin/env python3
import sys
import re

def parse_mulliken_charges(text):
  """
  Extract the LAST Mulliken charges block from Gaussian output.
  Returns list of (symbol, charge) tuples, e.g. ('N', 0.77734).
  """
  import re

  # Find ALL Mulliken sections
  matches = re.findall(
    r"Mulliken charges:(.*?)(Sum of Mulliken charges)",
    text,
    flags=re.S
  )

  if not matches:
    return []

  # Take ONLY the last Mulliken block
  block = matches[-1][0]

  charges = []
  line_re = re.compile(
    r"^\s*\d+\s+([A-Za-z]+)\s+([-+]?\d+\.\d+)",
    re.M
  )

  for sym, val in line_re.findall(block):
    charges.append((sym, float(val)))

  return charges
