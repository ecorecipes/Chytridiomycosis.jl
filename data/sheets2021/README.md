# Sheets et al. 2021: Bd viability (MTT) of five isolates at seven temperatures

- **File**: `mtt_od.csv`, the sheet `Ch1_MTT_OD` of the deposited workbook
  `Frontiers_Chytrid_2021_MTT_OD_Data.xlsx`, converted to CSV with values rounded to 6
  decimals (the workbook stores binary floating point, e.g. 0.11799999999999999 for 0.118).
  Missing readings are `NA`, as in the workbook. 5509 rows.
- **Source**: figshare, https://doi.org/10.6084/m9.figshare.14714841 ("Viability and Growth
  Dataset across Temperature Treatments"), file 28259625, retrieved 2026-11-19. The workbook,
  the zoospore-density workbook (https://doi.org/10.6084/m9.figshare.14714856, not bundled)
  and the authors' R Markdown (https://doi.org/10.6084/m9.figshare.14714865) are kept in
  `lab/papers/sheets2021/figshare/`.
- **Paper**: Sheets CN, Schmidt DR, Hurtado PJ, Byrne AQ, Rosenblum EB, Richards-Zawacki CL,
  Voyles J (2021). Thermal performance curves of multiple isolates of *Batrachochytrium
  dendrobatidis*, a lethal pathogen of amphibians. *Front Vet Sci* 8:687084.
  https://doi.org/10.3389/fvets.2021.687084
- **Licence**: CC BY 4.0. Attribution: Sheets et al. 2021, figshare.

## Columns (from the workbook's Metadata sheet)

| column | meaning |
|---|---|
| `Day` | day of the experiment after the zoospores were plated (0–43) |
| `Isolate` | `LA`, `TN`, `VT` (BdGPL1), `NM`, `OH` (BdGPL2) |
| `Temp` | incubation temperature (°C): 4, 12, 17, 21, 25, 26, 27 |
| `OD`, `ODNeg` | raw optical density of the live and heat-killed wells, no assay |
| `MTT`, `MTTNeg` | raw optical density (570 nm) of the live and heat-killed wells with the MTT assay |
| `MTTn` | `MTT − MTTNeg`: the normalised viability the logistic curves are fitted to |
| `ODn` | normalised OD |
| `Genotype` | `GPL1` or `GPL2` |
