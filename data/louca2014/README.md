# Louca et al. (2014): stored reproduction output

Louca S, Lampo M, Doebeli M (2014) Assessing host extinction risk following exposure to
*Batrachochytrium dendrobatidis*. *Proc R Soc B* 281:20132783.
<https://doi.org/10.1098/rspb.2013.2783>

**No data from the authors are bundled here.** The paper deposited no code or data. Its
electronic supplementary material (ESM, model details, Table S1, parameter estimates) was read
from the copy the first author's site appends to the article
(`https://pages.uoregon.edu/slouca/LoucaLab/SECTION_Publications/MODULE_Publications/Files/Louca2014Assessing.pdf`,
retrieved 2026-11-14; the journal's ESM link returned HTTP 403). That copy is © the authors /
The Royal Society and is kept outside the package (`lab/papers/louca2014_si/`). Every published
value the package uses is in `data/compendium/louca2014.toml`.

The *D. meridensis* infection series behind ESM Fig. S4 (Lampo 2013) are unpublished and were
not available.

## Files (our simulation output, MIT like the package)

- `reproduction.toml`: summary of our Monte Carlo reproduction (outcome fractions, risk effects,
  outcome probabilities binned by log10 Q_o and on a (β_A, μ_Z) grid) for the paper's four
  scenarios. Written by `scripts/louca2014_montecarlo.jl` (seed and design in the script).
- `montecarlo.csv`: one row per run: scenario, the scaled draws `u_*` (each in [0, 1] on its
  sampling scale; empty when not sampled in that scenario), `Qo`, `outcome`
  (`extinction`/`clearance`/`persistence`) and `t_end` (days).

These are outputs of our implementation, not published results.
