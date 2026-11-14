# Wilber et al. 2017: Bd transmission in a mesocosm density experiment

- **File**: `formatted_dd_data.csv`, 1 378 swab-to-swab transitions of adult *Rana muscosa*
  (*R. sierrae*) in 16 mesocosms (1, 4, 8 or 16 initially uninfected adults, each with five
  infected tadpoles), summer 2008.
- **Source**: Dryad, "Data from: Resistance, tolerance and environmental transmission
  dynamics determine host extinction risk in a load-dependent amphibian disease",
  https://doi.org/10.5061/dryad.p614h (version 1, published 2018-07-10), file
  `formatted_dd_data.csv` (140 860 bytes, MD5 b3c8668d7826a4cce4633b168fd11a11 as listed by
  the Dryad API). Dryad's download endpoint needed authentication (HTTP 401), so the copy in
  the authors' GitHub repository was used: https://github.com/mqwilber/host_extinction
  (master branch archive, file `data/formatted/formatted_dd_data.csv`; MD5
  b3c8668d7826a4cce4633b168fd11a11, identical to the Dryad listing; SHA-256
  6dd9c2ebf8225257c0bc050915b98cf1db7a16eb77f140fd01d850ea46675c9a). Retrieved 2026-11-14.
  The extracted archive is in `lab/papers/wilber2017_si/code/`.
- **Paper**: Wilber MQ, Knapp RA, Toothman M, Briggs CJ (2017) Resistance, tolerance and
  environmental transmission dynamics determine host extinction risk in a load-dependent
  amphibian disease. *Ecology Letters* 20:1169–1181. https://doi.org/10.1111/ele.12814
- **Licence**: CC0 1.0 (the Dryad record). The GitHub repository has no licence file. Only
  this Dryad-licensed file is bundled. The raw file `data/archival/DensExpt.csv` (all swabs,
  including tadpoles, needed to refit the latent zoospore-pool models of Table 2) is only in
  the repository and is **not** bundled.

## What was changed

Nothing.

## Fields

From the authors' `formatted_dd_data_metadata.txt` (Dryad README):

| Field | Meaning |
|:------|:--------|
| `size` | Bd load at the first swab, qPCR zoospore equivalents (ZE) per swab, natural scale; 0 = uninfected |
| `sizeNext` | load at the second swab (same units); empty when the frog died or was not swabbed |
| `swab_id` | id of the second swab |
| `time` | days between the swabs |
| `Day` | day of the experiment of the second swab (0 = start) |
| `individual` | frog id, `MDF<tank>_<frog>` |
| `density` | number of adults initially in the tank (1, 4, 8, 16) |
| `weight`, `svl`, `sex` | weight (g), snout–vent length (mm), sex (`F`, `M`, `M?`, or empty) |
| `ze_tadpoles` | total ZE on the tank's tadpoles at the first swab |
| `tank_num` | tank, 1–16 |
| `ze_frogs` | total ZE on the tank's frogs at the first swab |
| `num_infected` | number of infected frogs in the tank at the first swab |
| `cum_ze_frogs`, `cum_ze_tadpoles` | cumulative totals |

The transmission fits (Table 2, rows 1–3) use transitions with `Day < 35`, `3 < time < 9`,
no missing numeric field and `size == 0`: 333 transitions, 124 infections (the paper's
n = 333). Read it with `load_wilber2017()`.
