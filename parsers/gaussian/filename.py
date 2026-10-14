#!/usr/bin/env python3
import os
import re

METHOD_MAP = {
  "hf": "HF",
  "cam": "CAM-B3LYP",
  "camb3lyp": "CAM-B3LYP",
  "m062x": "M06-2X",
}

BASIS_MAP = {
  "631g": "6-31G",
  "631gdp": "6-31G(d,p)",
  "631++gdp": "6-31++G(d,p)",
  "6311g": "6-311G",
  "6311gdp": "6-311G(d,p)",
  "6311++gdp": "6-311++G(d,p)"
}

POPULATION_MAP = {
  "full": "full",
  "chelpg": "ChelpG",
  "hly": "HLY",
  "mk": "MK"
}

def parse_filename(fname):
  """Parse molecule-method-basis-population from a filename."""
  base = os.path.basename(fname)
  
  # strip extension
  name = re.sub(r"\.[A-Za-z0-9]+$", "", base)

  # split by - or _
  parts = re.split(r"[-_]", name)
  it = iter(parts)
  molecule = next(it)
  method = next(it, "")
  basis = next(it, "")
  population = next(it, "")

  method = METHOD_MAP.get(method, "")
  basis = BASIS_MAP.get(method, "")
  population = POPULATION_MAP.get(method, "")

  return {
    "filename": fname,
    "molecule": molecule,
    "method": method,
    "basis": basis,
    "population": population,
  }
