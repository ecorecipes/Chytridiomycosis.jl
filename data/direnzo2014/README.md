# DiRenzo et al. 2014: terminal infection intensity and zoospore output (Table 1)

- **File**: `table1.csv`, transcribed from Table 1 of the published article (one row per
  infected frog, in the printed order). No other individual-level values are published: the
  article has no supporting-information files and no data deposit, so the repeated swab and
  soak series used for the mixed models are not available.
- **Paper**: DiRenzo GV, Langhammer PF, Zamudio KR, Lips KR (2014). Fungal infection intensity
  and zoospore output of *Atelopus zeteki*, a potential acute chytrid supershedder. *PLOS ONE*
  9(3):e93356. https://doi.org/10.1371/journal.pone.0093356
- **Licence**: CC BY 4.0 (the article's licence). Attribution: DiRenzo et al. 2014, PLOS ONE.
- **Retrieved**: 2026-11-19, from the PLOS printable PDF.

## Columns

| column | meaning |
|---|---|
| `prior_exposure` | `Naive`, or `Previous` (exposed to the attenuated strain JEL 427-P39 23 weeks earlier) |
| `days_survived` | "Total days survived post-inoculation" |
| `intensity_at_death` | Bd infection intensity at death: zoospore genomic equivalents on one skin swab |
| `output_at_death` | zoospore output at death: zoospores released per minute in a 15-min soak |

## Known inconsistencies in the paper

- The abstract says all animals died between 22 and 33 days (mean 28.88, SE 1.58). The
  table's survival times run from 18 to 33 days, with mean 25.625 (SE 2.14). Figure 1 shows
  deaths at about 22, 25, 31, 32 and 33 days, which matches neither.
- The mean intensity at death (4,334,422) and its SE in the Results (1,156,576) are
  reproduced from the table; the abstract's SE (1,236,431) is not.
- The table's outputs give a mean of 23.29 (SE 21.07) and a maximum of 170.6; the text
  prints 23.55 (SE 22.78) and a range up to 172.61.
