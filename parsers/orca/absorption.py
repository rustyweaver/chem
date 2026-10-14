#!/usr/bin/env python3

import sys
import re

expr = (
  r"^----*\n"
  r"TD-DFT/TDA-EXCITATION SPECTRA\n"
  r"----*\n"
  r"\n"
  r"[^\n]*\n"
  r"[^\n]*\n"
  r"\n"
  r"----*\n"
  r"\s*ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\s*\n"
  r"----*\n"
  r"[^\n]*\n"
  r"[^\n]*\n"
  r"----*\n"
  r"(.*?\n)\n"
)

def parse_absorption_spectrum(text):
  """
  Parse the absorption spectrum lines
  """
  matches = re.findall(expr, text, flags=re.M|re.DOTALL)
  if not matches:
    result = []
  else:
    result = [
      {
        'from': cols[0],
        'to': cols[2],
        'energy': float(cols[3]),
        'wavelength': float(cols[5]),
        'fosc': float(cols[6]),
        'd2': float(cols[7]),
        'dx': float(cols[8]),
        'dy': float(cols[9]),
        'dz': float(cols[10])
      }
      for line in matches[-1].splitlines()
      if (cols := line.split()) and len(cols) >= 11
    ]
  return result
