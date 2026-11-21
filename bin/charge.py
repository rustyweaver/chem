#!/usr/bin/env python3

import sys
from spectrum_utils import get_charge

def fatal(s):
  print("ERROR: %s" % s, file=sys.stderr)
  sys.exit(1)

def main():
  if len(sys.argv) < 2:
    fatal("usage: charge.py <name>")
  print(get_charge(sys.argv[1]))

if __name__ == "__main__":
  main()
