# Midwife toads and a zoospore pool (Doddington et al. 2013)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the model of Mallorcan midwife toads (*Alytes
muletensis*) and free-living Bd zoospores of Doddington et al. (2013), *Ecology*
94:1795–1804 ([doi:10.1890/12-1270.1](https://doi.org/10.1890/12-1270.1)). The equations
(B.1–B.8) and parameters (Table B1) are in Appendix B; the stochastic model that produced the
paper's results is the authors' Berkeley Madonna code (Supplement 1). Both are on figshare
([doi:10.6084/m9.figshare.c.3306075](https://doi.org/10.6084/m9.figshare.c.3306075), CC0) and
were checked against the code here, together with the published version of the paper.

Time is in **years**, and zoospores are numbers in one pond. Tadpoles hatch from adults in the
breeding season, metamorphose in summer or overwinter, and juveniles mature the year after
metamorphosis. Zoospores infect tadpoles and overwintering tadpoles; infected tadpoles shed
zoospores; Bd kills only juveniles, at the extra rate `α`; infected juveniles mature into
uninfected adults.

```
dA/dt  = (1 − A/K) εJ (JX + JY) − μA A
dLX/dt = ψ A − (μL + εL + ν Z) LX           (L tadpoles, O overwintering tadpoles)
dOX/dt = − (μO + εO + ν Z) OX               (+ LX at time of year tM2)
dJX/dt = εL LX + εO OX − (μJ + εJ) JX
dLY/dt = ν Z LX − (μL + εL) LY
dOY/dt = ν Z OX − (μO + εO) OY              (+ LY at tM2)
dJY/dt = εL LY + εO OY − (μJ + α + εJ) JY
dZ/dt  = ρ (LY + OY) − μZ Z
```

`ψ`, `εL` and `εO` are nonzero only in their seasonal windows. Both problems use the code's
14 states, which split juveniles by age (this year's, and those able to mature) and by the
kind of tadpole they came from.

!!! note "Printed versus coded"
    The code matures juveniles at `εL` only during the tadpole metamorphosis window, and those
    from overwintering tadpoles at `εO` during theirs; Table B1 gives one constant
    `εJ = 12.3 yr⁻¹`. `maturation = :as_coded` (default) follows the code,
    `maturation = :as_printed` the table. The code's switch times (overwintering, juvenile
    ageing, Bd introduction) use strict inequalities that no time on its grid satisfies
    exactly; here each fires on the first step at or after its time.

## Results

- **Exposure trial (Table 1).** The group death rates and the isolate-specific excess rates
  with their 95% intervals are reproduced from the counts, except that the BdGPL rates come
  out 3.37 and 2.97 (published 3.38 and 2.98; Table B1 itself uses 2.97).
- **Fig. 3.** The stochastic model gives the pattern of Fig. 3: with high virulence
  (`α = 4.54`) host survival over 35 years falls to about 0.1 at the largest `ν` and `ρ` and
  is 1 at small values; with weak virulence (`α = 0.61`) it stays above about 0.75. The
  half-survival contour of Fig. 3a lies at about 2–3 times larger `ν` than the printed
  shading suggests. `scripts/doddington2013_fig3.jl` maps the full grid.
- **Most likely outcomes.** With `α = 4.54` Bd goes extinct after the epizootic, as
  published. With `α = 0.61` or `α = 0` at the Table B1 values of `ν` and `ρ`, Bd extinction
  is the most frequent outcome here, whereas the paper reports coexistence; coexistence is
  the modal outcome in about half of the Fig. 3b grid (the larger `ν` and `ρ`). These
  targets are recorded as discrepant.
- **Not implemented.** The SIR model with temperature-dependent recovery fitted in Fig. 4
  (the temperature series and counts are not deposited, and the recovery link is not
  printed), the fungal-load individual-based model (no code) and the state-space trend
  model of Appendix A.

```@docs
DoddingtonParameters
DoddingtonODEProblem
DoddingtonStochasticProblem
doddington_outcome
doddington_outcomes
doddington_grid
excess_death_rate
doddington_exposure_trial
```
