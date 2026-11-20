# Thermal performance of five Bd isolates (Sheets 2021)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the growth-curve analysis of Sheets et al. (2021),
*Front Vet Sci* 8:687084
([doi:10.3389/fvets.2021.687084](https://doi.org/10.3389/fvets.2021.687084)). Five BdGPL
isolates from the USA (LA, TN, VT: GPL1; NM, OH: GPL2) were grown at 4, 12, 17, 21, 25, 26 and
27 °C, and their viability followed by the MTT assay. For each isolate and temperature a
logistic curve is fitted to the normalised viability, by maximum likelihood with a standard
deviation that grows with the mean (the authors' R Markdown, figshare
[doi:10.6084/m9.figshare.14714865](https://doi.org/10.6084/m9.figshare.14714865)):

```
x(t) = K x0 / (x0 + (K − x0) exp(−r t))
MTTn_i ~ Normal(x(t_i), b0 + 1e-4 + b1 x(t_i))
```

The thermal response of an isolate is then `r` (or `K`) against temperature (Fig. 4, Table 2),
available as an [`IsolateTPC`](@ref) through [`sheets_tpc`](@ref). The data (figshare
[doi:10.6084/m9.figshare.14714841](https://doi.org/10.6084/m9.figshare.14714841), CC BY 4.0)
are bundled in `data/sheets2021`.

!!! warning "What is and is not reproduced"
    The refit reproduces 22 of the 35 fits of Table 2 to its four significant figures and four
    more to 2 %. Five fits (TN, OH and VT at 4 °C, TN at 26 °C, OH at 25 °C) differ from Table 2;
    the authors' own code, run on the deposited data, gives our values, so Table 2 came from
    other data or another run. K is not identified at LA 26 °C (printed 107300) or VT 27 °C, and
    VT 26 °C lies on a ridge in `(x0, r)`. Table 2's mean viabilities put the LA optimum at
    17 °C, not the 21 °C the Results state.

    The paper's Johnson–Lewin curves for `r` are in a Supplementary Material we could not
    retrieve; its Discussion calls them overfitted. They are not implemented.

```@docs
sheets_tpc
SheetsReading
sheets2021_path
load_sheets2021
SheetsGrowthModel
sheets_fits
```

[`maximum_likelihood`](@ref) and [`parameter_names`](@ref) also have methods for a
[`SheetsGrowthModel`](@ref).
