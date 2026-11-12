# Miller et al. 2025: *Litoria aurea* mesocosm state counts

- **Files**: `exp_data.csv` (the eight outdoor mesocosms of Waddle et al. 2024, as assigned
  to model states by Miller et al.) and `syn_data.csv` (the authors' synthetic data set,
  Supplement S2). Both unmodified.
- **Source**: Zenodo record "Chytrid Model v1.0", https://doi.org/10.5281/zenodo.15179344
  (published 2025-04-09), archive `clairemiller/chytrid_model-1.0.zip`, folder `data/`.
  Retrieved 2026-11-12. SHA-256:
  - `exp_data.csv` 0767da2cf2f9d12715457ee51798ad5574df07f3f09f0b0678efa8c54d240c42
  - `syn_data.csv` e9486bdfb7b6e7e0d4921aec9b11af685bf85dfee95560bd712342c15471bc51
- **Paper**: Miller CM, Lydeamore MJ, Flegg JA, Berger L, Skerratt LF, Waddle AW, Campbell PT
  (2025) Sunlight-heated refugia protect frogs from chytridiomycosis: a mathematical
  modelling study. arXiv:2503.06846 (v3, 5 May 2025; not peer reviewed as of 2026-11-12).
  https://doi.org/10.48550/arXiv.2503.06846
- **Underlying experiment**: Waddle AW et al. (2024) Hotspot shelters stimulate frog
  resistance to chytridiomycosis. *Nature* 631:344–349.
- **Licence**: CC BY 4.0 (the Zenodo record). Attribution: Miller et al. 2025.

## Fields

Long format, one row per mesocosm, week and state. The two files order the columns
differently, so read them by header.

| Field | Meaning |
|:------|:--------|
| `Mesocosm` | mesocosm id (exp: 2, 3, 5, 6, 8, 9, 11, 12; syn: 1–8) |
| `week` | weeks since the start: 0, 1, 2, 4, 6, 8, 10, 15 |
| `compartment` | `S`, `I1`, `I2`, `I3`, `SV`, `IV1`, `IV2`, `IV3` = the paper's S_U, I_U1–3, S_V, I_V1–3 |
| `shaded` | `TRUE` shade-cloth refugia; `FALSE` sunlight-heated ("unshaded") refugia |
| `N` | number of frogs in that state |

## Design (checked against the paper, §2.1.2 and §2.2.2)

Each mesocosm started with 20 frogs: 5 naive uninfected, 5 naive inoculated, 5 previously
infected and cleared ("vaccinated") uninfected, 5 vaccinated inoculated. Four mesocosms per
treatment (exp: shaded 3, 5, 9, 11; unshaded 2, 6, 8, 12). States were assigned from each
frog's infection intensity at the current and previous swab (paper Table 2), with the
high-load threshold 10⁵ (naive) or 10⁴ (vaccinated). Frogs not found are absent from the
counts, so totals fall (and occasionally rise again, when a missed frog is re-found; 3.4% of
observations). At week 0 three mesocosms (2, 3, 6) record inoculated vaccinated frogs in
`IV2` rather than `IV1`; the model always starts from `IV1 = 5`.
