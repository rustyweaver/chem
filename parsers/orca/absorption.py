#!/usr/bin/env python3
import sys
import re


"""
----------------------------------------------------------------------------------------------------
                     ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS    
----------------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ   
                      (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)  
----------------------------------------------------------------------------------------------------
  0-1A  ->  1-1A    3.977164   32078.0   311.7   0.302218378   3.10163   1.73796  -0.10381   0.26521
  0-1A  ->  2-1A    4.002773   32284.5   309.7   0.065529072   0.66821   0.80938  -0.03821   0.10800
  ...

"""


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
    r"(?P<body>(?:[^\n]*\n)*?)"   # ← CAPTURE all lines non-greedily
    r"\n\n"                       # ← ends with two blank lines
)

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
        'fosc': float(cols[6])
      }
      for line in matches[-1].splitlines()
      if (cols := line.split())
    ]
  return result
