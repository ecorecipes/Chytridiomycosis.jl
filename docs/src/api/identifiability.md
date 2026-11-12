# Swabs, eDNA and identifiability

```@meta
CurrentModule = Chytridiomycosis
```

!!! note "New model (N5)"
    This observation model and the identifiability results are new work in this package,
    not published results. The assay scales used in examples are illustrative assumptions.

An observation model that links the Briggs et al. (2010) load model
([doi:10.1073/pnas.0912886107](https://doi.org/10.1073/pnas.0912886107)) to skin-swab qPCR
(per frog, or pooled) and to eDNA in pond water, and tools to ask which Briggs parameters each
data stream can recover.

In the linear (sub-lethal) phase, the mean sporangial load per frog `m`, the zoospore pool `Z`
and the excess variance of load among frogs `w = Var − m` obey a linear system
([`BriggsMomentProblem`](@ref)); `k = m²/w` is the aggregation of loads (negative-binomial
shape, as in [`LoadBDI`](@ref)). A swab tube of `n` pooled swabs is negative binomial with
shape `n k` and mean `n ρ m`; an eDNA sample is lognormal around `κ Z/V`.

What each design identifies (structural, by hand and confirmed by
[`sensitivity_rank`](@ref)):

| Data | Identified |
|---|---|
| swab means, one density | trace and determinant of the linear system only |
| swab means (or eDNA alone), two or more densities | `γ/V`, `μ`, `r0 = ηνf − σ`, `νη(1−f)` |
| + eDNA, both assays on absolute scales | `ν` and `η(1−f)` |
| + aggregation `k` from individual swabs | `ηf`, so every parameter |

Without an absolute swab scale (copies per sporangium), the symmetry `ρ → cρ`, `ν → ν/c`
remains. At one density, `γ` and `μ` are never separated.
[`expected_information`](@ref) measures how well a design does in practice; fits use
[`maximum_likelihood`](@ref) and [`profile_likelihood`](@ref), which have methods for
[`SurveyModel`](@ref), and [`parameter_names`](@ref).

## Latent moments

```@docs
BriggsMomentProblem
```

## Design, data and likelihood

```@docs
SurveyPond
SurveySample
SurveyObservation
simulate_survey
SurveyModel
```

## Identifiability

```@docs
identifiable_combinations
sensitivity_rank
expected_information
observed_information
```
