#!/usr/bin/env python3
# split_goat_xyz.py
# (Original code preserved above...)

# --- Extract and print the Final Ensemble Info block from ORCA output ---
import sys
from pathlib import Path
import re

if __name__ == "__main__":
    if len(sys.argv) > 1:
        out_file = Path(sys.argv[1])
        if out_file.exists():
            content = out_file.read_text()
            pattern = re.compile(
                r"#\s*Final ensemble info\s*#.*?(?:\n\s*-+\n)(.*?)(?=\n\s*\n|\Z)",
                re.DOTALL
            )
            match = pattern.search(content)
            if match:
                print("=== Final Ensemble Info ===")
                print(match.group(0).strip())
            else:
                print("No 'Final ensemble info' section found.")
        else:
            print(f"File not found: {out_file}")
    else:
        print("Usage: python split_goat_xyz.py goat.out")
