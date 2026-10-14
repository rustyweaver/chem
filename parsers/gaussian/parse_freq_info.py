#!/usr/bin/env python3

import sys
import re
import numpy as np

# Regex for Gaussian scientific notation (D) and markers
derivative_expr = r"Dipole derivative wrt mode\s+\d+:\s+(\S+)\s+(\S+)\s+(\S+)"
inertia_expr = r"Principal axes and moments of inertia in atomic units:"

def parse_derivatives_array(text):
  """
  Return Nx3 array of dipole derivatives from Gaussian output text.
  """
  matches = re.findall(derivative_expr, text)
  if not matches:
    return np.array([])
  
  # Replace 'D' with 'e' and convert to float
  data = [[float(val.replace('D', 'e')) for val in row] for row in matches]
  return np.array(data)


def print_derivatives(data_array):
  if data_array.size == 0:
    print("derivatives = np.array([])")
    return
    
  arr_str = ",\n  ".join([str(list(r)) for r in data_array])
  print(f"derivatives = np.array([\n  {arr_str}\n])")



def parse_last_inertia_axes(text):
  """
  Extracts the 3x3 eigenvector matrix from the final inertia block.
  Skips header lines and parses columns 1, 2, and 3.
  """
  matches = list(re.finditer(inertia_expr, text))
  if not matches:
    return np.array([])
  
  # Skip the 'Principal' line and the following two lines (Eigenvalues and X label)
  start_idx = matches[-1].end()
  lines = text[start_idx:].split('\n')[3:6]
  
  axes = []
  for line in lines:
    parts = line.split()
    # Parse columns 1, 2, and 3 (indices 1, 2, and 3)
    if len(parts) >= 4:
      axes.append([float(x) for x in parts[1:4]])
  return np.array(axes)



def print_eigenvectors(axes_array):
  if axes_array.size == 0: return
  arr_str = ",\n  ".join([str(list(r)) for r in axes_array])
  print(f"eigenvectors = np.array([\n  {arr_str}\n])")



def parse_rotational_constants(text):
  """
  Return a map of rotational constants from the last MHz block in Gaussian output.
  Uses only readable keys: 'equilibrium', 'ground_vibrational', 'ground_distorted'
  """
  blocks = re.findall(r'Rotational Constants \(in MHz\).*?Ce=.*?C00=.*?C0=.*?(?=\n\s*-+|$)', text, re.S)
  if not blocks:
    return {}

  last_block = blocks[-1]
  
  # Map internal labels to readable keys
  mapping = {'e': 'equilibrium', '00': 'ground_vibrational', '0': 'ground_distorted'}
  res = {v: {'A': 0.0, 'B': 0.0, 'C': 0.0} for v in mapping.values()}

  # Capture Axis (A,B,C), Label (e,00,0), and Value
  pattern = r'([ABC])(e|00|0)=\s+([\d.]+)'
  matches = re.findall(pattern, last_block)
  
  for axis, label, value in matches:
    res[mapping[label]][axis] = float(value)

  return res

def print_rotational_constants(mhz_data_map):
  """
  Print the rotational constants map using only readable keys.
  """
  
  def format(dmap):
    order = ['equilibrium', 'ground_vibrational', 'ground_distorted']
    lines = []
    for key in order:
      vals = dmap[key]
      line = f"  '{key}': {{'A': {vals['A']}, 'B': {vals['B']}, 'C': {vals['C']}}}"
      lines.append(line)
    return ",\n".join(lines)

  formatted = format(mhz_data_map)
  print("# Rotational constants (MHz)")
  print(f"rotational_constants_mhz = {{\n{formatted}\n}}")

  conv = lambda x: (x * 1e6) / 29979245800.0
  cm_data_map = {k: {sub_k: conv(sub_v) for sub_k, sub_v in v.items()} for k, v in mhz_data_map.items()}
  formatted = format(cm_data_map)
  print("# Rotational constants (cm-1)")
  print(f"rotational_constants_cm = {{\n{formatted}\n}}")


def parse_frequencies(text):
  header_pattern = r"Fundamental Bands\s*\n\s*-+\s*\n\s*Mode\(n\)\s+E\(harm\)\s+E\(anharm\)\s+I\(harm\)\s+I\(anharm\)"
  
  assert re.search(header_pattern, text), (
    'Parsing error look for the Anharmonic Infrared Spectroscopy header.')

  # Extract the lines following the header
  data_section = re.split(header_pattern, text)[1].strip()
  results = []

  for line in data_section.splitlines():
    match = re.search(r"^\s*\d+\(\d+\)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)", line)
    if match:
      results.append({
        "eharm": float(match.group(1)),
        "eanharm": float(match.group(2)),
        "iharm": float(match.group(3)),
        "ianharm": float(match.group(4))
      })
    elif not line.strip():
      break
    else:
      assert False, 'Parsing error'

  results.sort(key=lambda x: x['eharm'])
  return results


def print_frequencies(data_list):
  lines = []
  for entry in data_list:
    line = f"  {{'eharm': {entry['eharm']}, 'eanharm': {entry['eanharm']}, 'iharm': {entry['iharm']}, 'ianharm': {entry['ianharm']}}}"
    lines.append(line)

  formatted = ",\n".join(lines)
  print("# Anharmonic Infrared Spectroscopy (E in cm^-1, I in km.mol^-1)")
  print(f"anharmonic_frequencies = [\n{formatted}\n]")


if __name__ == "__main__":
  text = sys.stdin.read()
  derivatives = parse_derivatives_array(text)
  eigenvectors = parse_last_inertia_axes(text)
  rotational_constants = parse_rotational_constants(text)
  frequencies = parse_frequencies(text)

  assert len(derivatives) > 0
  assert len(eigenvectors) > 0
  assert len(rotational_constants) > 0
  assert len(frequencies) > 0

  print_derivatives(derivatives)
  print_eigenvectors(eigenvectors)
  print("projection = derivatives @ eigenvectors")
  print_rotational_constants(rotational_constants)
  print_frequencies(frequencies)
