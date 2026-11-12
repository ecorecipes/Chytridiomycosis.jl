# Thermal performance

```@meta
CurrentModule = Chytridiomycosis
```

Temperature-dependent rates for the mechanistic models. A curve maps temperature (°C) to a
rate. [`IsolateTPC`](@ref) and [`HostTPC`](@ref) attach a curve to a pathogen isolate or a
host trait, with its provenance. [`SeasonalTemperature`](@ref) and
[`temperature_dependent`](@ref) turn a curve into a time-varying rate.

Published thermal data are in the [compendium](compendium.md):

- **Cardinal temperatures.** Limits and optima from Stevenson et al. 2013 and Sheets et
  al. 2021 (Bd), and Martel et al. 2013 (Bsal). There are no fitted curves.
- **Kelly et al. 2024.** Fitted growth-rate polynomials for 13 Bsal isolates
  (`IsolateTPC(:kelly2024, "Eifel2019")`).
- **Woodhams et al. 2008.** The Bd life-history rates at several temperatures.
- **Blooi et al. 2015.** Heat clearance of Bsal (25 °C for 10 days).

A curve built from cardinal temperatures, such as a [`Briere`](@ref) through Martel's
values, rests on the choice of family and is not a published curve.

## Curves

```@docs
ThermalPerformanceCurve
Briere
CardinalTPC
GaussianQuadratic
SharpeSchoolfield
Ratkowsky
PolynomialTPC
InterpolatedTPC
thermal_optimum
thermal_data
```

## Isolates, hosts and forcing

```@docs
IsolateTPC
HostTPC
SeasonalTemperature
temperature_dependent
```

## Bd life history in culture (Woodhams et al. 2008)

```@docs
WoodhamsParameters
intrinsic_growth_rate
```
