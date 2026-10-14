#!/usr/bin/env python3
import sys
import re

def parse_dipole_moment(text):
  """
  Extract the LAST Dipole moment, in debeye, from gaussian output that
  looks like this:

  Electric dipole moment (dipole orientation):
  (Debye = 10**-18 statcoulomb cm , SI units = C m)
                   (au)            (Debye)         (10**-30 SI)
    Tot        0.667558D+00      0.169676D+01      0.565979D+01
  """

  # Find all dipole moment sections and take the last. Debye is the
  # third column.
  matches = re.findall(
    r"Electric dipole moment.*\sTot\s+\S+\s+(\S+)",
    text,
    flags=re.S
  )
  if not matches:
    return []
  dptext = matches[-1]

  # Convert fortran D to python E
  dptext = dptext.replace('D', 'E')
  return float(dptext)
