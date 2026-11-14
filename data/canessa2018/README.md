# Canessa et al. 2018: archived parameter draws for the Bsal IPM

- **File**: `parameter_draws.csv`, the 1,000 random parameter sets used for the integral
  projection model simulations (one row per simulation).
- **Source**: Dryad, "Data from: Decision making for mitigating wildlife diseases: from
  theory to practice for an emerging fungal pathogen of amphibians",
  https://doi.org/10.5061/dryad.0jh5p (published 2018-12-18), file
  `Canessa_JAPPL_2017_00978_SI2.xls` (298,496 bytes, MD5 29341e442ab6f46801108c74973bce89,
  "Full set of random values used for integral projection model simulations (parameter
  distributions)"). The same file is `lab/papers/canessa2018_si4.xls` (identical MD5). The
  Dryad record lists no other files. Metadata retrieved 2026-11-13.
- **Paper**: Canessa S, Bozzuto C, Campbell Grant EH, et al. (2018) Decision-making for
  mitigating wildlife diseases: from theory to practice for an emerging fungal pathogen of
  amphibians. *J Appl Ecol* 55:1987-1996. https://doi.org/10.1111/1365-2664.13089
- **Licence**: CC0 1.0 (the Dryad record).

## What was changed

The single sheet (`beta0`) was read with xlrd 2.0.1 and written out with the shortest
round-trip decimal for every IEEE double, so the values are exactly those in the `.xls`
(a LibreOffice conversion was checked and rejected because it rounds to 15 digits). A `draw`
column (row number, 1-1000) was added. The sheet's three trailing columns, which hold the
authors' descriptions of the columns, were dropped; their text is quoted below.

## Fields

| Field | Meaning (authors' description in quotes) |
|:------|:--------|
| `draw` | row number, 1-1000 (added) |
| `S_0` … `S_10000` | "Load-dependent survival for a three-day interval, for infection loads of 0, 10, 100, 1000, 10000 GE (values drawn from aggregate distribution of elicited values)" |
| `T_0` … `T_10000` | "Load-dependent trasmission (probability of transmitting infection to a susceptible individual) for infection loads of 0, 10, 100, 1000, 10000 GE" (same source); `T_0` is 0 in every row |
| `L_max` | maximum load on an infected host (GE) ("Maximum load and growth rate of Bsal on infected individual host (values drawn from aggregate distribution of elicited values)") |
| `L_growth` | zoospores produced by one zoospore in five days (same description; the paper, SI lines 109-111) |
| `beta0` | "Contact rate of salamanders derived from Stegen et al. 2017" (ha per 3-day step) |
| `S0` | "Density of salamanders (n/ha) derived from Stegen et al. 2017"; every value is an integer count divided by 0.475 ha (475 m transect × 10 m) |

Loads are Bsal genomic equivalents (GE, qPCR). Read with `load_canessa2018()`.
