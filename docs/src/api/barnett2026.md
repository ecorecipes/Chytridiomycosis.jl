# Imperfect prophylaxis ABM (Barnett et al. 2026, preprint)

```@meta
CurrentModule = Chytridiomycosis
```

!!! warning "Preprint"
    This model comes from a **bioRxiv preprint that has not been peer reviewed**: Barnett et
    al. (2026), "From protection to amplification: imperfect chytridiomycosis prophylaxis
    increases infections in wild amphibians",
    [doi:10.64898/2026.05.15.725113](https://doi.org/10.64898/2026.05.15.725113) (v2,
    August 2026). The agent-based model is a port of the authors' NetLogo code
    ([kmbarnett100/protection_to_amplification](https://github.com/kmbarnett100/protection_to_amplification)),
    checked against all their released runs.

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
(`barnett_vaccine(...; adverse = true)`). Every value (Table S1, Methods and the NetLogo code) is in
`data/compendium/barnett2026.toml`.

## The code and the text

The port follows the code where it differs from the preprint's text. The source header of
`src/models/barnett2026.jl` lists every difference. The ones that matter most are these:

- **Infection.** Zoospores on a patch are shared among all hosts in the world, not only those
  on the patch. Only the share that falls on hosts of that patch infects them. A vaccinated
  host's lower establishment probability reduces the number of establishing zoospores for
  everybody, but not that host's own share.
- **Clearance** is a rate: daily survival of a zoosporangium is ``e^{-0.2} = 0.82``.
- **Stalled infections.** A host whose zoosporangia number 0 stops maturing the stages it
  carries until it is infected again.
- **Population.** Every perimeter patch starts with 500 tadpoles, one infected. Ponds have
  about 136 patches, and 28–49 of them are on the perimeter.
- **Vaccine response** varies by ±5%, not ±10%, and not at all in the two-way runs.
- **Statistics.** Intensity is the median load of infected metamorphs on the last day, and
  changes are relative to the median control run. Runs with fewer than two infected metamorphs
  are left out.

## What is reproduced

`scripts/barnett2026_abm.jl` runs the package with the parameters of each of the authors'
29,325 runs, eight times over, and compares the distributions (`data/barnett2026/`). NetLogo's
random numbers cannot be reproduced, so the match is in distribution.

| Scenarios (mean over runs: authors / package) | Survivors | Median load | Zoospore density |
| --- | --- | --- | --- |
| Control | 128.8 / 128.4 | 261.6 / 252.6 | 399,160 / 392,056 |
| Anti-infection | 208.3 / 206.8 | 85.7 / 86.9 | 220,203 / 219,712 |
| Anti-transmission | 193.6 / 193.6 | 102.9 / 102.7 | 246,809 / 246,563 |
| Anti-growth | 187.8 / 187.2 | 125.9 / 126.1 | 288,713 / 287,356 |
| Tolerance | 147.3 / 149.3 | 315.0 / 314.2 | 436,701 / 432,909 |
| Tolerance × anti-infection | 220.3 / 219.6 | 71.8 / 71.7 | 193,627 / 192,298 |
| Tolerance × anti-transmission | 210.8 / 210.4 | 85.4 / 85.9 | 218,437 / 220,594 |
| Tolerance × anti-growth | 206.8 / 206.3 | 115.6 / 114.8 | 271,688 / 268,308 |
| Adverse: susceptibility | 40.1 / 41.3 | 300.2 / 304.2 | 520,100 / 522,678 |
| Adverse: shedding | 42.9 / 48.5 | 312.0 / 305.1 | 558,665 / 548,248 |
| Adverse: clearance | 53.9 / 53.1 | 307.2 / 303.6 | 478,661 / 482,744 |
| Adverse: tolerance | 92.5 / 96.5 | 238.1 / 229.6 | 321,888 / 315,239 |

Of 85 comparisons of means (up to nine outputs per file), 77 have |z| < 2 and one has
|z| > 3: the survivors of the increased-shedding runs (42.9 against 48.5, z = 3.6), recorded as
discrepant. That gap does not grow with efficacy, and is not explained. The adverse-tolerance
file shows smaller differences in the same direction (|z| up to 2.9). The other ten files agree
within sampling error.

Every statement of the preprint about its scenario figures that the package checks has the
published sign, by a clear margin. This includes three that an earlier version of this page
recorded as discrepant, when only the text was available.

!!! note "Increased susceptibility"
    The preprint reports that infection intensity rises when vaccination increases
    susceptibility (Fig. 5a). That holds where metamorphs survive (about +25% at 50% efficacy
    and coverage). At high efficacy and coverage almost none do (none of the authors' 75 runs at
    +100% and coverage of at least 80%), so that corner of the figure is extrapolated.

## Field trial and challenge experiment

[`BarnettBACIModel`](@ref) refits the authors' four mixed models (zero-inflated negative
binomial for intensity, binomial for prevalence; pond and year as random effects), with the
Laplace likelihood that glmmTMB uses. Log-likelihoods, estimates, standard errors and marginal
means agree with glmmTMB to four or five digits. That includes the headline interaction
(P = 0.00135) and the 88% and 23% declines. The preprint's text prints the intervals of the two
control-pond prevalences swapped. The ratio it gives as 6.42 is 6.44 before rounding.

!!! note "Data not bundled"
    The repository has no licence, so its files are not bundled. [`barnett2026_path`](@ref)
    finds them in the lab workspace or in `ENV["CHYTRIDIOMYCOSIS_BARNETT2026"]`.

## API

```@docs
BarnettParameters
BarnettVaccine
barnett_vaccine
BarnettABMProblem
barnett_summary
barnett_relative_change
barnett2026_path
load_barnett2026
BarnettBACIModel
maximum_likelihood(::BarnettBACIModel)
```
