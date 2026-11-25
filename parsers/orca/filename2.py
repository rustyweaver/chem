#!/usr/bin/env python3

import os
import re
import sys

MOLECULES = [
  "nbi-327",
  "nbi-526",
  "daota",
  "atota",
  "tota",
  "tata",
  "aadota",
  "adota",
]

METHODS = [
  "hf",
  "cam-b3lyp",
  "b3lyp",
  "xtb",
  "pbe0",
]

BASIS_SETS = [
  "def2-tzvpd",
  "631pgs",
]

SOLVENTS = [
  "dcm",
  "water",
  "mecn",
  "toluene",
  "ch3cn",
  "vacuum",
]

def check_keywords(keywords):
  n = len(keywords)
  for i in range(n):
    for j in range(i + 1, n):
      assert keywords[i] not in keywords[j], "%s in %s" % (keywords[i], keywords[j])

def get_keyword(s, keywords):
  for keyword in keywords:
    if keyword in s:
      return keyword
  return ""

check_keywords(MOLECULES)
check_keywords(METHODS)
check_keywords(BASIS_SETS)
check_keywords(SOLVENTS)

def parse_filename(filename):
  """Parse molecule-method-basis-solvent from a filename."""
  base = os.path.basename(filename).lower()
  
  # strip extension
  name = re.sub(r"\.[A-Za-z0-9]+$", "", base)

  attributes = {
    'filename' : filename,
    'molecule' : get_keyword(name, MOLECULES),
    'solvent' : get_keyword(name, SOLVENTS),
    'method' : get_keyword(name, METHODS),
    'basis' : get_keyword(name, BASIS_SETS),
  }
  return attributes


if __name__ == "__main__":
  for arg in sys.argv[1:]:
    print('%s: %s' % (arg, parse_filename(arg)))
