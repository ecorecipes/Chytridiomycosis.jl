# Turner et al. 2025: survival series and start populations

- **Source**: Zenodo, "Data and code for 'Timing of environmental watering influences
  extinction risk for chytrid-infected amphibian populations'",
  https://doi.org/10.5281/zenodo.10921748 (`TurnerBdModel.zip`, 200 MB). Retrieved
  2026-11-19 and kept outside the package in `lab/papers/turner2025/zenodo/`.
- **Paper**: Turner A, Heard G, Mathwin R, Bradshaw CJA, Hall A, Wassens S (2025).
  Cool-season environmental water delivery increases extinction risk for chytrid-infected
  amphibians. *Conserv Sci Pract* 7:e70022. https://doi.org/10.1111/csp2.70022
- **Licence**: CC BY 4.0. Attribution: Turner et al., Zenodo.

## Files (extracts, not the original files)

The archive's inputs are nine 40 MB arrays `referenceData/<Month>WateringArray_<pattern>.csv`
(10,000 replicates × 52 years × 6: wetness, adult survival mean and SD, metamorph summer
survival mean and SD, an empty column) and nine `Terminal_structure_<Month>_<pattern>.csv`
start populations. Two small extracts are bundled:

| file | content |
|---|---|
| `survival_series.csv` | for each month (`august`, `september`, `october`) and dry-year pattern (`random` = RandomDryYears, `max2dry` = 2DryYears, `max1dry` = 1DryYears), the 52 yearly values of adult survival mean and SD and summer survival mean and SD from the array's LAST replicate (row 10,000), which is what the authors' `pop.func` uses for every replicate (it indexes with `iter`, not `rpt`). |
| `start_populations.csv` | the first 1000 non-empty start populations (eggs, females aged 1–5) of each scenario's `Terminal_structure` file. Files with 10,002 columns carry an extra row-number column (values 1..6) after the index; it is dropped. |

The extraction was done with R (`read.csv`, `array(…, dim = c(10000, 52, 6))`, as the authors
read the arrays). Wetness sequences are not bundled: `turner_watering_sequence` draws them
by the rule the arrays follow, or `load_turner2025_survival_array` reads a full array from the
extracted archive.
