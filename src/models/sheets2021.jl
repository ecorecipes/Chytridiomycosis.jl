# Sheets, Schmidt, Hurtado, Byrne, Rosenblum, Richards-Zawacki & Voyles (2021) Front Vet Sci
# 8:687084, doi:10.3389/fvets.2021.687084. Five BdGPL isolates (LA, TN, VT: GPL1; NM, OH: GPL2)
# grown at 4, 12, 17, 21, 25, 26 and 27 °C. Checked against the full text (PMC8258153) and the
# authors' R Markdown (figshare doi:10.6084/m9.figshare.14714865).
#
# For each isolate and temperature, a logistic curve is fitted to the normalised MTT
# viability series (Supplementary R code):
#
#   x(t) = K x0 / (x0 + (K − x0) exp(−r t))     (dX/dt = r X (1 − X/K), X(0) = x0)
#   MTTn_i ~ Normal(x(t_i), b0 + 1e-4 + b1 x(t_i))
#
# by maximum likelihood (multi-start Nelder–Mead in optim()). Time in days. The thermal
# response of each isolate is then r (and K) against temperature (Fig. 4, Table 2). The paper
# also fitted Johnson–Lewin curves to r, in a Supplementary Material we could not retrieve
# (Data_Sheet_1.PDF); its Discussion calls those fits overfitted, and they are not implemented.

const SHEETS_ISOLATES = ("LA", "TN", "VT", "NM", "OH")
const SHEETS_TEMPERATURES = (4.0, 12.0, 17.0, 21.0, 25.0, 26.0, 27.0)

"""
    sheets_tpc(isolate; trait = :r, extrapolation = :zero) -> IsolateTPC

The thermal response of one Bd isolate of Sheets et al. (2021) (`"LA"`, `"TN"`, `"VT"`,
`"NM"` or `"OH"`), as an [`IsolateTPC`](@ref) that interpolates linearly
([`InterpolatedTPC`](@ref)) between the seven temperatures of Table 2. `trait` is

- `:r`: the logistic growth rate of normalised MTT viability (day⁻¹);
- `:K`: its carrying capacity (OD);
- `:MTT_mean`: the mean maximum viability (OD at 570 nm);
- `:zoospore_density`: the mean maximum zoospore density (Table 4; no value at 25 °C).

The values are the paper's (compendium `:sheets2021`, sets `"<isolate>_<T>C"`). The authors
show (Fig. 4) that `r` stays high at 25–27 °C while `K` collapses, so `r` alone is a poor
summary of growth at high temperature; and K at LA 26 °C (printed 107300) is not identified.
"""
function sheets_tpc(isolate::AbstractString; trait::Symbol = :r, extrapolation::Symbol = :zero)
    isolate in SHEETS_ISOLATES ||
        throw(ArgumentError("isolate must be one of $(SHEETS_ISOLATES); got \"$isolate\""))
    trait in (:r, :K, :MTT_mean, :zoospore_density) ||
        throw(ArgumentError("trait must be :r, :K, :MTT_mean or :zoospore_density"))
    Ts, vs = thermal_data(:sheets2021, trait; set_prefix = isolate * "_")
    rec = first(r for r in parameter_records(:sheets2021) if startswith(r.set, isolate * "_"))
    label = trait === :r ? "logistic growth rate r of normalised MTT viability (day⁻¹)" :
            trait === :K ? "logistic carrying capacity K of normalised MTT viability (OD)" :
            trait === :MTT_mean ? "mean maximum viability, MTT optical density at 570 nm" :
            "mean maximum zoospore density"
    return IsolateTPC("Bd", String(isolate), get(rec.system, "lineage", ""), label,
        InterpolatedTPC(Ts, vs; extrapolation), :sheets2021)
end
