# Bsal in fire salamanders (Schmidt et al. 2017)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the susceptible–latent–infectious model of Bsal in fire
salamanders of Schmidt et al. (2017), *R Soc Open Sci* 4:160801
([doi:10.1098/rsos.160801](https://doi.org/10.1098/rsos.160801)). It is checked against the
paper. Parameters are read from the [compendium](compendium.md):
`SchmidtParameters(:schmidt2017)`.

Four values are not reproduced (see `data/validation/targets.toml`): S_eq and K_H below,
and the shares of the German habitat infected after 8 years by scattered colonies (Results;
ESM fig. S3), where eq. (2.5) divided by the habitat area gives 8.97 and 0.128 against the
paper's ≈9% and ≈0.1% (see [`scattered_colony_area`](@ref)).

- **S_eq.** The paper gives 1.0027. Eq. (2.3) gives 1.0030. The paper's value equals
  `(d + e)(d + d_I)/(βe)`, the large-K limit of eq. (2.3), so the authors probably reported
  that approximation.
- **K_H.** The paper gives about 9.5. The Jacobian gives 10.21. No rate, unit or model
  variant tested gives 9.5. The equilibrium is barely stable below 10.21: the real part of
  the eigenvalues is −0.0011 yr⁻¹ at K = 9.5.

```@docs
SchmidtParameters
SchmidtSLIProblem
invasion_threshold
removal_fraction
basic_reproduction_number
endemic_equilibrium
hopf_threshold
wave_speed
scattered_colony_area
population_crash
```

## Spatial spread and long-distance jumps

The reaction–diffusion version (§2.2) is discretised by finite volumes and lowered to an
`ODEProblem`; its front speed converges to [`wave_speed`](@ref). The scattered-colony model
behind eq. (2.5) is simulated as a jump process with colonies placed in a square habitat, so
that overlaps can be counted.

!!! note "Modifications of the published model"
    The planar reaction–diffusion model is the paper's. The radial geometry
    (`SchmidtPDEProblem(p; geometry = :radial)`) and the explicit colony simulation
    ([`SchmidtColonyProblem`](@ref)), which counts overlaps and the habitat edge that eq. (2.5)
    ignores, are modifications made in this package.

```@docs
SchmidtPDEProblem
front_position
front_speed
SchmidtColonyProblem
```
