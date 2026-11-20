# Brannelly et al. 2018: *Acris crepitans* capture–recapture and Bd swab data

- **Paper**: Brannelly LA, Chatfield MWH, Sonn J, Robak M, Richards-Zawacki CL (2018). Fungal
  infection has sublethal effects in a lowland subtropical amphibian population. *BMC
  Ecology* 18:34. https://doi.org/10.1186/s12898-018-0189-5 (open access, CC BY 4.0;
  PMC6137908).
- **Licence of the data**: the article states that "the Creative Commons Public Domain
  Dedication waiver (http://creativecommons.org/publicdomain/zero/1.0/) applies to the data
  made available in this article, unless otherwise stated". **CC0 1.0.**
- **How obtained**: the three Additional files, downloaded on 2026-11-19 from the Europe PMC
  supplementary-files endpoint for PMC6137908 (full copy in
  `lab/papers/brannelly2018/supp/`). The only change is the line endings: the originals use
  carriage returns (old Mac style) and were converted to `\n`. Cell contents are unchanged.

| file here | original | contents |
|---|---|---|
| `additional_file1_captures.csv` | `12898_2018_189_MOESM1_ESM.csv` | Additional file 1: one row per capture event of *A. crepitans* (699 rows, the 699 capture events of the Results, Nov 2010–Jul 2012): date, animal ID, recapture flag, SVL, mass, body and microhabitat temperature, ZE, log10 ZE, infection status, disease category |
| `additional_file2_other_species.csv` | `12898_2018_189_MOESM2_ESM.csv` | Additional file 2: swabs of the other species (species code, month, year, log ZE, status) |
| `additional_file3_capture_histories.csv` | `12898_2018_189_MOESM3_ESM.csv` | Additional file 3: capture histories of the 315 marked animals over the 19 survey days of Nov 2010–Nov 2011, with SVL, mass, temperatures and body condition |

## Coding of Additional file 3

The header says: 0 not captured, 1 captured and Bd-negative, 2 captured and Bd-infected.
The file also has 9 entries coded `3`, not described. Matching them with Additional file 1
(animal and date) shows they are infected captures in disease category 2 (log10 ZE above
2.455), while `2` is category 1. With `3` counted as infected, the file gives 378 capture
events, 315 negative and 63 positive, the totals of Table 1 (whose "Positive" and
"Negative" labels are swapped against the text, which says 63 captures were Bd-positive).
Table 1 also gives 266 unique individuals; the file has 315 animals, of which 50 were caught
on more than one survey day (Table 1: 50 recaptured).

The columns to the right of the 19 dates contain a few spreadsheet remnants (`#REF!`, a
stray `10`); the loader ignores them.
