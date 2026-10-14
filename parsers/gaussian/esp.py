#!/usr/bin/env python3
import sys
import re

def parse_esp_charges(text):
  """
  Extract the LAST ESP charge block from Gaussian output.
  Returns list of (symbol, charge) tuples, e.g. ('N', -1.115833).
  """
  import re

  # Find ALL ESP sections
  matches = re.findall(
    r"ESP charges:(.*?)(Sum of ESP charges)",
    text,
    flags=re.S
  )

  if not matches:
    return []

  # Take ONLY the last block
  block = matches[-1][0]

  charges = []
  line_re = re.compile(
    r"^\s*\d+\s+([A-Za-z]+)\s+([-+]?\d+\.\d+)",
    re.M
  )

  for sym, val in line_re.findall(block):
    charges.append((sym, float(val)))

  return charges
