# Bozzuto & Canessa 2019: R0(t) curves digitised from Fig. 2

- **File**: `fig2_R0t.csv`, the grey R0(t) curves of Fig. 2a (hibernating population) and
  Fig. 2b (aestivating population), digitised here. These are not data from the authors.
- **Source**: Bozzuto C, Canessa S (2019). Impact of seasonal cycles on host-pathogen
  dynamics and disease mitigation for *Batrachochytrium salamandrivorans*. *Glob Ecol
  Conserv* 17:e00551. https://doi.org/10.1016/j.gecco.2019.e00551. The figure is the
  embedded 300 dpi raster on p. 5 of the publisher's PDF (`lab/papers/bozzuto2019.pdf`).
  Digitised 2026-11-14.
- **Licence**: the article is CC BY-NC-ND 4.0. The values are measurements taken from its
  figure, with attribution: Bozzuto & Canessa 2019, Fig. 2.

## Method

The raster was read in greyscale. The y scale comes from the centres of the tick labels:
in panel (a), 0 at pixel row 525.0 and 21.6 px per unit; in panel (b), 0 at row 1153.5 and
87.5 px per unit. The x axis runs from month 1 to the end of month 12, pixel columns 184 to
637, mapped to t = 0–1 yr. For each column, the value is the darkness-weighted centre of the
longest run of grey pixels (grey level 150–238), after masking the black Floquet-R0 and R0
lines. One pixel is 0.05 units in (a) and 0.011 units in (b).

## Columns

| column | meaning |
|---|---|
| `panel` | `a` (hibernation) or `b` (aestivation) |
| `t` | time (decimal years, 0 = start of January) |
| `R0t` | R0(t), read from the grey curve |

Used in `test/bozzuto2019.jl` to check which aestivation input explains the 1.6% offset of
the aestivation results (see `bozzuto2019_R0_aestivation` in `data/validation/targets.toml`).
