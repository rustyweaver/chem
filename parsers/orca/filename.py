#!/usr/bin/env python3

import os
import re
import sys

ALIAS_MAP = {
  r"(^|-)(cam-b3lyp-d4)(-|$)" : r"\g<1>camb3lyp\g<3>",
  r"(^|-)(6-31\+\+g--)(-|$)" : r"\g<1>631++g**\g<3>",
  r"^(restart-|)dha-conf-" : r"dha",
  r"^(restart-|)vhf-conf-" : r"vhf",
  r"(^|-)aug-cc-pvdz(-|$)" : r"\g<1>augccpvdz\g<2>",
  r"(^|-)aug-pc-0(-|$)" : r"\g<1>augpc0\g<2>",
  r"(^|-)aug-pc-1(-|$)" : r"\g<1>augpc1\g<2>",
  r"(^|-)aug-pc-2(-|$)" : r"\g<1>augpc2\g<2>",
  r"(^|-)aug-pc-3(-|$)" : r"\g<1>augpc3\g<2>",
  r"(^|-)pc-0(-|$)" : r"\g<1>pc0\g<2>",
  r"(^|-)pc-1(-|$)" : r"\g<1>pc1\g<2>",
  r"(^|-)pc-2(-|$)" : r"\g<1>pc2\g<2>",
  r"(^|-)pc-3(-|$)" : r"\g<1>pc3\g<2>",
  r"(^|-)m06-2x(-|$)" : r"\g<1>m062x\g<2>",
  r"neb-xtb-ts" : r"nebxtbts",
  r"neb-dft-ts" : r"nebdftts",
  r"neb-dft-with-guess-ts" : r"nebdftwithguessts",
  r"def2-tzvp" : r"def2tzvp",
}

METHOD_MAP = {
  "hf": "HF",
  "cam": "CAM-B3LYP",
  "camb3lyp": "CAM-B3LYP",
  "m062x": "M06-2X",
}

BASIS_MAP = {
  "631g": "6-31G",
  "631gdp": "6-31G(d,p)",
  "631++g": "6-31++G",
  "631++gdp": "6-31++G**",
  "631++g**": "6-31++G**",
  "6311g": "6-311G",
  "6311gdp": "6-311G(d,p)",
  "6311++gdp": "6-311++G(d,p)",
  "augccpvdz": "aug-cc-pVDZ",
  "augpc0": "aug-PC-0",
  "augpc1": "aug-PC-1",
  "augpc2": "aug-PC-2",
  "augpc3": "aug-PC-3",
  "pc0": "PC-0",
  "pc1": "PC-1",
  "pc2": "PC-2",
  "pc3": "PC-3",
  "def2tzvp": "def2-TZVP",
}

POPULATION_MAP = {
  "full": "full",
  "chelpg": "ChelpG",
  "hly": "HLY",
  "mk": "MK"
}

SOLVENT_MAP = {
  "dcm": "DCM",
  "mecn": "MeCN",
  "toluene": "Toluene",
}

MOLECULE_MAP = {
  "nebxtbts" : "neb-xtb-ts",
  "nebdftts" : "neb-dft-ts",
  "nebdftwithguessts" : "neb-dft-with-guess-ts",
}

def parse_filename(filename):
  """Parse molecule-method-basis-solvent from a filename."""
  base = os.path.basename(filename).lower()
  
  # strip extension
  name = re.sub(r"\.[A-Za-z0-9]+$", "", base)

  for pattern, replacement in ALIAS_MAP.items():
    name = re.sub(pattern, replacement, name)


  # split by - or _
  parts = re.split(r"[-_]", name)

  attributes = {
    'filename' : filename,
    'molecule' : parts[0],
    'solvent' : ''
  }

  maps = {
    'molecule' : MOLECULE_MAP,
    'method' : METHOD_MAP,
    'basis' : BASIS_MAP,
    'pop' : POPULATION_MAP,
    'solvent' : SOLVENT_MAP,
  }

  for part in parts:
    matches = 0
    for (attribute, map) in maps.items():
      if part in map:
        attributes[attribute] = map[part]
        matches += 1
    assert matches <= 1, 'Ambiguous attribute "%s" in filename %s' % (part, filename)

  return attributes


if __name__ == "__main__":
  for arg in sys.argv[1:]:
    print('%s: %s' % (arg, parse_filename(arg)))
