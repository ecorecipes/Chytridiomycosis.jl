# Imperfect prophylaxis ABM (Barnett et al. 2026, preprint)

```@meta
CurrentModule = Chytridiomycosis
```

!!! warning "Preprint"
    This model comes from a **bioRxiv preprint that has not been peer reviewed**: Barnett et
    al. (2026), "From protection to amplification: imperfect chytridiomycosis prophylaxis
    increases infections in wild amphibians",
    [doi:10.64898/2026.05.15.725113](https://doi.org/10.64898/2026.05.15.725113) (v2,
    August 2026). The authors' NetLogo code was not released, so this implementation was
    checked against the text only.

The preprint couples a field trial with an agent-based model. In the field trial, Bd
metabolites (a prophylaxis that acts like an imperfect vaccine) were added to whole ponds in
California. The agent-based model asks how such an intervention changes Bd in a cohort of
Pacific chorus frogs (*Pseudacris regilla*). It follows tadpoles over a 90-day season,
through metamorphosis, on a pond made of three kinds of patch: perimeter (shallow), deep and
terrestrial. Each host carries a zoosporangium load. Zoospores pass between hosts through
the perimeter patches, and hosts also reinfect themselves. Natural mortality is daily and
stage-specific. A metamorph dies of Bd when its load reaches `smax = 562`.

A vaccinated tadpole has one or more of its traits changed:

| Mode | Trait changed | Protective direction |
|---|---|---|
| anti-infection resistance (`:infection`) | establishment probability `est` | decrease |
| anti-transmission resistance (`:transmission`) | zoospores shed per zoosporangium | decrease |
| anti-growth resistance (`:growth`) | daily clearance of zoosporangia | increase |
| tolerance (`:tolerance`) | lethal load `smax` | increase |

The preprint's "adverse reaction" scenarios reverse each direction
(`barnett_vaccine(...; adverse = true)`). Every published value (Table S1 and Methods) is in
`data/compendium/barnett2026.toml`.

## Choices where the preprint is silent

The source header of `src/models/barnett2026.jl` lists every choice. The ones that matter
most are these:

- **Population and pond size.** The preprint gives neither. The package uses 1000
  tadpoles, so that 0.2% initially infected is 2 tadpoles, and a ring of 60 perimeter
  patches.
- **Transmission.** The preprint describes a multinomial draw driven by an exposure
  parameter (0.25) and zoospore mortality (2 day⁻¹). The package reads these as competing
  daily hazards.
- **Vaccine efficacy.** A change of `e` is applied as the factor `(1 ± e)^imm`, which is
  the preprint's exponential modulation. Here `imm` is uniform on 1 ± 10%.
- **The intensity statistic.** The preprint does not define it. The default
  (`metric = :metamorph`) is the season mean of the mean load of infected live metamorphs.
  `:season` uses all live hosts, and `:final` uses the survivors on the last day.

## What is reproduced

The targets in `data/validation/targets.toml` test the signs of the preprint's percent
changes against an unvaccinated control.

!!! note "Replicate runs"
    The preprint ran 25 replicates per scenario. With 25, the run-to-run spread of a percent
    change (standard deviations of 5–37 percentage points across seeds) is as large as several
    of the effects, so their signs depend on the seed and on the Julia version's random
    streams. The targets therefore use 200 replicates, and 1,000 for the two smallest effects.

The model reproduces the preprint's headline results:

- Resistance (anti-infection, anti-transmission or anti-growth) lowers infection intensity.
  Anti-infection resistance also increases the final population. The anti-growth effect is
  small (−7% to −14% with 1,000 replicates).
- **Tolerance raises infection intensity.** This is the amplification that the preprint
  offers as one explanation of its field result.
- **Increased susceptibility, an adverse reaction, lowers the population.** This is the
  second explanation. Its effect on intensity is only weak with the default statistic
  (about +2%, with a sign that is not robust across seeds).
- Reduced tolerance lowers intensity.

Three statements are **discrepant**:

- **Resistance does not counteract tolerance** in the two-way scenarios (Fig. 4a, Fig. S7):
  with tolerance, strong anti-infection resistance gives about +40% and anti-growth
  resistance about +63%, where the preprint reports decreases.
- The **control's median prevalence** is 0, not 100%.

Self-reinfection (17.8 × 0.1 × 0.25 ≈ 0.45 per zoosporangium per day) outgrows clearance
even when clearance is doubled. As a result, infected tadpoles reach the lethal load by
metamorphosis. With the season-mean statistic, the two-way results have the preprint's sign,
but tolerance alone no longer does. The field-trial numbers are `pending`, because the field
data are not deposited.

## API

```@docs
BarnettParameters
BarnettVaccine
barnett_vaccine
BarnettABMProblem
barnett_summary
barnett_relative_change
```
