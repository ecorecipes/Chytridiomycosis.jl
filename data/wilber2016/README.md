# Wilber et al. 2016: Bd loads of adult frogs in a temperature experiment

- **File**: `converted_temp_data.csv`, 589 swab-to-swab transitions of 20 adult frogs.
- **Source**: Dryad, "Data from: Integral Projection Models for host-parasite systems with an
  application to amphibian chytrid fungus", https://doi.org/10.5061/dryad.06sd9 (version 1,
  published 2017-03-11), file `converted_temp_data.csv` (18 103 bytes, MD5
  9d29dfec21b278bc3293bccef8df2b66 as listed by the Dryad API). Dryad's download endpoint
  needed authentication, so the byte-identical copy in the authors' GitHub repository was
  used: https://github.com/mqwilber/ipm_for_parasites (commit 240b5e67, file
  `data/formatted/converted_temp_data.csv`; MD5 9d29dfec21b278bc3293bccef8df2b66, SHA-256
  8d0eddc981a87e2aee9bcfbd6715dabbb1277b515393d5dcb3e1cc50021cfcfe). Retrieved 2026-11-13.
  The repository archive, extracted, is in `lab/papers/wilber2016_si/repo/`.
- **Paper**: Wilber MQ, Langwig KE, Kilpatrick AM, McCallum HI, Briggs CJ (2016) Integral
  projection models for host–parasite systems with an application to amphibian chytrid
  fungus. *Methods Ecol Evol* 7:1182–1194. https://doi.org/10.1111/2041-210X.12561
- **Licence**: CC0 1.0 (the Dryad record). The GitHub repository has no licence file; only the
  Dryad-licensed data file is bundled here, not the code or its results.

## What was changed

Nothing. The file is the authors' `converted_temp_data.csv`, produced by their
`convert_temperature_data.py` from the raw `TempExp.csv` (also in the repository, not bundled).

## Fields

| Field | Meaning |
|:------|:--------|
| `size` | Bd load at the first swab, qPCR zoospore equivalents (ZE) per swab, natural scale; 0 = uninfected |
| `sizeNext` | load at the next swab (same units); empty when the frog died between the swabs |
| `time` | days between the swabs: 3, or 8 for the first interval after exposure, or 120 for a few gaps |
| `individual` | frog id, `F<temperature>_<k>` |
| `temp` | temperature treatment, °C: 4, 12, 20 or 26 |
| `surv` | 1 if the frog was alive at the second swab, 0 if it died |
| `stage` | always `adult` |

The species is given as *Rana sierrae* in the raw file; the paper calls the complex *R.
muscosa*. Frogs were exposed to about 10⁶ zoospores and swabbed every 3 days from day 8 to
day 119. The paper's analyses use only 3-day intervals and exclude 26 °C (five frogs), which
the text does not mention ("20 adult *R. muscosa* … at three different temperatures … 5 frogs
per temperature"). Read it with `load_wilber2016()`.
