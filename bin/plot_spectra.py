#!/usr/bin/env python3

# Plot the spectra from csv files.

import argparse
import csv
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from scipy.signal import find_peaks
from spectrum_utils import parse_keywords
from spectrum_utils import keyword_match
from spectrum_utils import get_display_name
from spectrum_utils import reformat
import sys

parser = argparse.ArgumentParser(description="Plot absorption spectra.")
parser.add_argument("--outfile", required=True, help="Output PNG file")
parser.add_argument("--infile", nargs="+", required=True, help="Input CSV files")
parser.add_argument("--images", nargs="+", default=[], help="Molecular PNG images")
parser.add_argument("--min_wavelen", type=int, default=285)
parser.add_argument("--max_wavelen", type=int, default=600)
parser.add_argument("--molecule", default="")
parser.add_argument("--method", default="")
parser.add_argument("--basis", default="")
parser.add_argument("--solvent", default="")
parser.add_argument("--sample", default="")
parser.add_argument("--source", default="")
parser.add_argument("--spectrum_type", default="")
parser.add_argument("--title", default="")
parser.add_argument("--legend", default="solvent")

args = parser.parse_args()

fields = ["molecule", "method", "basis", "solvent", "sample", "source", "spectrum_type"]
filters = []

for field in fields:
    for filter in getattr(args, field).split(','):
        if filter:
            filters.append((field, filter))

legend_fields = args.legend.split(",")

if any(field not in fields for field in legend_fields):
    parser.error("--legend must contain molecule, method, basis, or solvent")

spectra = []

for infile in args.infile:
    print('Reading', infile)

    with open(infile, newline="") as f:
        reader = csv.DictReader(f)

        if not reader.fieldnames or reader.fieldnames[0] != "wavelength":
            parser.error(f"{infile}: first column must be 'wavelength'")

        rows = list(reader)

    wl = np.array([float(row["wavelength"]) for row in rows])

    # Take the source and spectrum_type only from the filename
    source = parse_keywords(infile)['source']
    spectrum_type = parse_keywords(infile)['spectrum_type']

    for column in reader.fieldnames[1:]:
        values = parse_keywords(column)
        values['source'] = source
        values['spectrum_type'] = spectrum_type
        if not keyword_match(values, filters):
            continue

        intensity = np.array([float(row[column]) if row[column] else np.nan for row in rows])
        label = " / ".join(get_display_name(values[field]) for field in legend_fields if values[field])

        spectra.append((wl, intensity, label, source, spectrum_type))

if not spectra:
    sys.exit("ERROR: No matching spectra found, check command-line arguments.")


# The max value for normalization is performed separately for
# experimental vs. calculated, because these use different scales.

source_max = {}

for wl, intensity, label, source, spectrum_type in spectra:
    mask = (wl >= args.min_wavelen) & (wl <= args.max_wavelen)

    if np.any(np.isfinite(intensity[mask])):
        maximum = np.nanmax(intensity[mask])
        source_max[(source, spectrum_type)] = max(source_max.get((source, spectrum_type), 0), maximum)

if not source_max or any(value <= 0 for value in source_max.values()):
    sys.exit("No positive intensities for one or more sources")


fig = plt.figure(figsize=(10, 3))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1])

# fig = plt.figure(figsize=(9, 4))
# gs = fig.add_gridspec(1, 2, width_ratios=[2, 1])
ax = fig.add_subplot(gs[0, 0])

right = gs[0, 1].subgridspec(2, 1, height_ratios=[1, 1])
ax_table = fig.add_subplot(right[0])
ax_table.axis("off")

ax_images = fig.add_subplot(right[1])
ax_images.axis("off")

rows = []

for wl, intensity, label, source, spectrum_type in spectra:
    if label == "":
        label = "Other"

    normalized = intensity / source_max[(source, spectrum_type)]
    ax.plot(wl, normalized, label=get_display_name(label))

    if spectrum_type != 'absorption':
        continue

    mask = np.isfinite(normalized) & (wl >= args.min_wavelen) & (wl <= args.max_wavelen)

    x = wl[mask]
    y = normalized[mask]

    peaks, _ = find_peaks(y, prominence=0.05)
    peak_wavelengths = x[peaks]
    peak_wavelengths = np.sort(peak_wavelengths)[::-1]

    if len(peak_wavelengths) >= 2:
        w1 = peak_wavelengths[0]
        w2 = peak_wavelengths[1]

        e1 = 1239.84 / w1
        e2 = 1239.84 / w2

        rows.append([
            label,
            f"{w1:.1f}",
            f"{w2:.1f}",
            f"{w1-w2:.1f}",
            f"{e2-e1:.3f}"
        ])

if len(rows) > 0:
    table = ax_table.table(
        cellText=rows,
        colLabels=["Spectrum", "Peak 1 (nm)", "Peak 2 (nm)", "Δλ (nm)", "ΔE (eV)"],
        colWidths=[0.29, 0.13, 0.13, 0.13, 0.13],
        cellLoc="center",
        loc="upper left",
    )

    table.auto_set_font_size(False)
    table.set_fontsize(6)
    table.scale(1, 1.2)

if args.title:
    ax.set_title(reformat(args.title), fontsize=7)

ax.set_xlim(args.min_wavelen, args.max_wavelen)
ax.set_ylim(0.0, 1.05)
ax.set_xlabel("Wavelength (nm)", fontsize=7)
ax.set_ylabel("Normalized Absorbance / Intensity (au)", fontsize=7)
ax.legend(frameon=False, fontsize=7)
ax.xaxis.set_minor_locator(MultipleLocator(10))

fig.tight_layout()

x = 0.0

for filename in args.images:
    image = plt.imread(filename)
    zoom = 0.16
    artist = OffsetImage(image, zoom=zoom)
    ax_images.add_artist(AnnotationBbox(
        artist,
        (x, 0.5),
        xycoords="axes fraction",
        box_alignment=(0, 0.5),
        frameon=False,
        pad=0
    ))

    x += image.shape[1] * zoom / (fig.dpi * ax_images.get_position().width * fig.get_figwidth())

print('Writing', args.outfile)

dirname = os.path.dirname(args.outfile)

if dirname:
    os.makedirs(dirname, exist_ok=True)

fig.savefig(args.outfile, dpi=300)
plt.close(fig)
