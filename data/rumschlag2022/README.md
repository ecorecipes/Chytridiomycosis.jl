# Rumschlag et al. 2022: zoospore survival and tadpole infection assays

- **Files**: the deposited CSV files, unmodified except for their names. The deposit's
  names are swapped relative to their contents.
  - `survival.csv` is the deposit's `TadpoleData.csv`: the zoospore survival assay, 199
    units.
  - `infection.csv` is the deposit's `BdSurvival.csv`: the tadpole infection assay, 197
    tadpoles.
- **Source**: figshare, https://doi.org/10.6084/m9.figshare.16814701.v1 (posted
  2021-10-14), files 31090369 and 31090372. Retrieved 2026-11-12. The original files and the
  metadata spreadsheet are kept in `lab/data/rumschlag2022/`.
- **Paper**: Rumschlag SL, Roth SA, McMahon TA, Rohr JR, Civitello DJ (2022). Variability in
  environmental persistence but not per capita transmission rates of the amphibian chytrid
  fungus leads to differences in host infection prevalence. *J Anim Ecol* 91:170–181.
  https://doi.org/10.1111/1365-2656.13612
- **Licence**: CC BY 4.0. Attribution: Rumschlag et al. 2021, figshare.

## Water codes (deposit metadata)

`ASW` = artificial spring water, `A` = Trout Creek, `B` = Flatwoods, `C` = Green Swamp.

## Columns

**`survival.csv`**

| column | meaning |
|---|---|
| `time.point` | nominal sampling time: 0, 6, 16, 24 or 48 h |
| `time.hrs` | recorded time of the count (h) |
| `live.cell.count`, `dead.cell.count` | haemocytometer counts (Trypan blue) |
| `live.cells.ul`, `dead.cells.ul` | cells per µl (count × 3.33325) |

**`infection.csv`**

| column | meaning |
|---|---|
| `time.trt` | lag before the tadpole was added: 0, 3, 6, 12 or 24 h |
| `body.mass`, `mouth.mass` | masses (g) |
| `load.TM` | qPCR load |
| `infect.TM` | infected (1) or not (0) |

- Each tadpole was exposed for one day in 100 ml.
- Host: *Osteopilus septentrionalis*. Bd isolate: SRS 810. Temperature: 21 °C.
- Read by `load_rumschlag2022()`.
