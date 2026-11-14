# Zoospore persistence (Rumschlag et al. 2022)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the zoospore persistence model of Rumschlag et al.
(2022), *J Anim Ecol* 91:170–181
([doi:10.1111/1365-2656.13612](https://doi.org/10.1111/1365-2656.13612)). Its equations and
values were checked against the **published version**. The bioRxiv preprint
([doi:10.1101/2021.03.17.435818](https://doi.org/10.1101/2021.03.17.435818)) prints the same
numbers.

The refit uses the authors' deposited assay data (figshare,
[doi:10.6084/m9.figshare.16814701](https://doi.org/10.6084/m9.figshare.16814701), CC BY 4.0),
bundled in `data/rumschlag2022`. It reproduces the paper's estimates:

- the lifespans 14.4, 13.3, 11.4 and 9.4 h;
- the persistence time of 21.6 h;
- their "± 95% CI" values;
- the likelihood ratio χ² of 7.581;
- 13 of the 15 AICc values in Table 1.

The other two are the models in which both `m` and `β` have two estimates (ASW and ponds).
Ours are 3.3 units lower. Table 1's own values show that the original fits of these two
models fell short of the maximum. Each has a lower implied log-likelihood, by 1.32 and 1.37
units, than the same model with a single `β`, which is nested in it. No other nested pair in
Table 1 is out of order.

## Model

Bd zoospores were followed in four waters at 21 °C: artificial spring water (ASW) and water
from three Florida ponds (Trout Creek, Flatwoods, Green Swamp). Cuban treefrog tadpoles were
the hosts. Time is in days, and zoospores are counted per µl.

```
dZL/dt = −m ZL           live zoospores
dZD/dt =  m ZL − d ZD    dead zoospores
dS/dt  = −β S ZL         susceptible tadpoles
dI/dt  =  β S ZL         infected tadpoles
```

The paper reports the lifespan `24/m` and the persistence time `24/d` in **hours**.

!!! warning "Two corrections to the paper"
    - **Eq. 6.** As printed, eq. 6 has its bracket reversed and is negative for every
      `t > 0`. [`pool_state`](@ref) uses the correct solution.
    - **The scale of β.** The paper gives no units for β. The data identify only the
      product `β ZL(0)`. The printed `β = 7.1 × 10⁻⁶` equals that product divided by `10⁵`,
      so it reproduces the paper's Fig. 4 only with `ZL0 = 1e5` in
      [`infection_probability`](@ref). That is an inferred convention.

```@docs
RumschlagParameters
ZoosporePersistenceProblem
zoospore_lifespan
persistence_time
infection_probability
```

## A live/dead zoospore pool for other models

The live and dead compartments form a reusable component. A qPCR or eDNA count sees DNA from
live and dead zoospores (`ZL + ZD`), but only live zoospores infect (`ZL`). Other models can
add this pool to their own right-hand side through [`pool_derivatives`](@ref).

```@docs
ZoosporePool
pool_derivatives
pool_state
qpcr_zoospores
infectious_zoospores
viable_fraction
```

## Data and refit

!!! note "Reproduction"
    The refit reproduces a published analysis; it is not a new model. Survival times are
    the recorded times of the counts, because Table 1 is not reproduced with the nominal
    times. One initial concentration `ZL0` is fitted and shared by both assays. That gives
    the paper's degrees of freedom.

```@docs
RumschlagCount
RumschlagTadpole
rumschlag2022_path
load_rumschlag2022
RumschlagModel
rumschlag_parameters
rumschlag_model_selection
rumschlag_covariance
rumschlag_summary
```

[`maximum_likelihood`](@ref) and [`parameter_names`](@ref) also have methods for a
[`RumschlagModel`](@ref).
