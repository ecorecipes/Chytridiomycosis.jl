# Pradel recruitment and growth (Hardy 2023)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the capture–recapture analysis of Hardy et al. (2023),
*J Appl Ecol* 60:121–131
([doi:10.1111/1365-2664.14310](https://doi.org/10.1111/1365-2664.14310)): a robust-design
Pradel model of adult male boreal toads (*Anaxyrus boreas boreas*) at five Colorado sites,
1998–2020, with Bd arrival as a before/after covariate. The data are the authors' Dryad
deposit ([doi:10.5061/dryad.fttdz08wv](https://doi.org/10.5061/dryad.fttdz08wv), CC0),
bundled in `data/hardy2023`.

!!! note "Version read"
    Equations and numbers were checked against the version of record and its Supporting
    Information (Appendices S1.1–S1.6). They were first read from the reprint in Chapter 1 of
    B. M. Hardy's PhD dissertation (Colorado State University, 2023), which prints the same
    numbers; only its table and figure numbers differ.

The fit reproduces:

- the capture counts of Table 1 and the 1,594 marked toads;
- 26 of the 27 focal models of Appendices S1.3–S1.5 to 0.01 in −2 log L (and their K and AICc);
- the post-Bd survival (0), recruitment (0.25, SE 0.02) and growth rate (0.25, 95% CI
  0.21–0.30);
- the range of capture probabilities (0.03–0.81) and the pre-Bd survival at Morgan's Gulch
  (0.83, SE 0.01).

Nine targets are discrepant. Our fit of φ(t^Bd) is 9.28 units of −2 log L better than MARK's.
The top model's AICc is 26649.35 against 26649.10. The printed −2 log L of Table 3
(10503.97) does not go with its own AICc and K. The second model of Table 3 is 48.45 AICc
units behind, not 56.81. The pre-Bd λ range (0.69–1.40 against 0.65–1.46), the elevation
coefficients (0.27 and −0.11 against 0.28 and −0.12) and the survival at Delilah Pond (0.48
against 0.47) differ by less than one printed SE. See the `hardy2023_*` targets.

Muths, Scherer & Pilliod (2011, *J Appl Ecol*,
[doi:10.1111/j.1365-2664.2011.02005.x](https://doi.org/10.1111/j.1365-2664.2011.02005.x))
used the same reverse-time model on two boreal toad populations and suggested compensatory
recruitment after Bd; Hardy et al. test that idea. Its data are not deposited, and it is not
implemented here.

## The Pradel temporal-symmetry model

Pradel (1996) models a population forward in time with apparent survival `φ` and backward
in time with seniority `γ`. In the `f` parameterization the per-capita recruitment is `f`,
the realized growth rate is `λ = φ + f`, and `γ_{t+1} = φ_t / λ_t`. In the robust design,
each primary period (year) has several secondary occasions with capture probability `p`
(here `p = c`), and `p* = 1 − ∏(1 − p)`. The likelihood is conditional on the number of
animals ever caught, and abundance is derived as `N = n / p*`. The layer is general (any
number of groups, any design matrices, standard or robust design). It reproduces Program
MARK's −2 log L exactly.

```@docs
PradelGroup
pradel_loglikelihood
pradel_capture_probability
pradel_abundance
PradelProblem
PradelFit
pradel_derived
pradel_se
pradel_coefficient
simulate_pradel
```

## Boreal toads in Colorado (Hardy et al. 2023)

Two choices were made where the paper is silent. Both are checked by the printed −2 log L.

- An interval `t → t + 1` takes the Bd status of year `t`.
- The climate covariates for recruitment are those of three years before the recruits are
  counted.

The drought and winter-onset columns of two of the deposited covariate files carry swapped
headers. The loader assigns them by value.

```@docs
hardy2023_path
Hardy2023Site
load_hardy2023
hardy2023_captures
hardy2023_problem
hardy2023_fit
hardy2023_model_set
hardy2023_summary
```
