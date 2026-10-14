# flatten.py
"""
Convert parsed charge data (nested lists of atoms/charges)
into long-format atom-level rows suitable for Pandas.
"""

CANONICAL_POPULATION = {
  "full":    "Mulliken",
  "ChelpG":  "ESP/ChelpG",
  "MK":      "ESP/MK",
  "HLY":     "ESP/HLY",
}


def flatten_by_charge(record):
  out = []

  for model in ("mulliken", "esp"):
    charges = record.get(model)
    if not charges:
      continue

    assert "population" in record, "bad record - missing population"
    population = record["population"]
    assert population in CANONICAL_POPULATION, "bad population %s" % population
    canonical_population = CANONICAL_POPULATION[population]

    for i, (atom, charge) in enumerate(charges, start=1):
      out.append({
        "molecule": record["molecule"],
        "method": record["method"],
        "basis": record["basis"],
        "population": canonical_population,
        "atom_index": i,
        "atom": atom,
        "charge": charge
      })

  return out

def flatten_all_by_charge(records):
  """
  Apply flatten_by_charge to a list of parse records.
  Returns a flat unified list of atom-level rows.
  """
  flat = []
  for rec in records:
    flat.extend(flatten_by_charge(rec))
  return flat
