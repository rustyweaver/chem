#!/usr/bin/env python3

from __future__ import annotations
import re

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
  "cam-b3lyp-d4",
  "b3lyp",
  "xtb"
]

BASIS_SETS = [
  "def2-tzvpd",
]

SOLVENTS = [
  "dcm",
  "water",
  "h2o",
  "mecn",
  "toluene",
  "ch3cn",
  "vacuum",
  "etoh",
]

SAMPLES = [
  "s1",
  "s2",
  "s3",
  "s4",
  "s5",
  "s6",
]

SOURCES = [
  "experimental",
  "calculated",
]

SOURCES = [
  "experimental",
  "calculated",
]

SPECTRUM_TYPES = [
  "absorption",
  "emission",
]

DISPLAY_NAMES = {
  "nbi-327" : "NBI-327",
  "nbi-526" : "NBI-526",
  "daota" : "DAOTA+",
  "atota": "ATOTA+",
  "tota": "TOTA+",
  "tata": "TATA+",
  "adota": "ADOTA+",
  "aadota": "AADOTA+",
  "ch3cn" : "CH3CN",
  "water" : "H2O",
  "h2o" : "H2O",
  "dcm" : "DCM",
  "vacuum" : "Vacuum",
  "etoh" : "EtOH",
  "b3lyp" : "CAM-B3LYP",
  "def2-tzvpd" : "def2-TZVPD",
  "s1" : "Sample 1",
  "s2" : "Sample 2",
  "s3" : "Sample 3",
  "s4" : "Sample 4",
  "s5" : "Sample 5",
  "s6" : "Sample 6",
  "calculated" : "Calculated",
  "experimental" : "Experimental",
}

def check_keywords(keywords):
  n = len(keywords)
  for i in range(n):
    for j in range(i + 1, n):
      assert keywords[i] not in keywords[j], "%s in %s" % (keywords[i], keywords[j])

def get_keyword(s, keywords):

  def normalize(keyword):
    return "h2o" if keyword == "water" else keyword

  for keyword in keywords:
    if re.search(r"(^|-|/)" + re.escape(keyword) + r"(-|/|$)", s):
      return normalize(keyword)
  return ""

check_keywords(MOLECULES)
check_keywords(METHODS)
check_keywords(BASIS_SETS)
check_keywords(SOLVENTS)
check_keywords(SAMPLES)
check_keywords(SOURCES)
check_keywords(SPECTRUM_TYPES)

def reformat(s):
  return (s.replace('CH3CN', '$\mathrm{CH}_3\mathrm{CN}$')
          .replace('H2O', '$\mathrm{H}_2\mathrm{O}$'))


# Returns the name to display.
def get_display_name(keyword):
  display_name = DISPLAY_NAMES.get(keyword, keyword)
  return reformat(display_name)


# Parse molecule-method-basis-solvent from a string.
def parse_keywords(name):
  keywords = {
    'molecule' : get_keyword(name, MOLECULES),
    'solvent' : get_keyword(name, SOLVENTS),
    'method' : get_keyword(name, METHODS),
    'basis' : get_keyword(name, BASIS_SETS),
    'sample' : get_keyword(name, SAMPLES),
    'source' : get_keyword(name, SOURCES),
    'spectrum_type' : get_keyword(name, SPECTRUM_TYPES),
  }
  return keywords


def keyword_match(values, filters):
  for key, filter in filters:
    if not values.get(key):
      continue
    if filter.startswith('!'):
      if values[key] == filter[1:]:
        return False
    else:
      if not values[key] == filter:
        return False
  return True


# Apply subscripting to selected chemical names in strings used in the
# figures.
def reformat(s):
  return (s.replace('CH3CN', '$\mathrm{CH}_3\mathrm{CN}$')
          .replace('H2O', '$\mathrm{H}_2\mathrm{O}$'))


# Find map keys that contain the pattern.
def find_matching_keys(map, pattern: str | list[str]):
  patterns = [pattern] if isinstance(pattern, str) else pattern
  return [k for k in map if any(re.search(p, k) for p in patterns)]


def get_solvent_of(soln, solvents):
  for solv in solvents:
    match_string = ' in ' + solv
    if match_string in soln:
      return solv
  assert False, 'Unable to determine solvent for "' + soln + '"' 
