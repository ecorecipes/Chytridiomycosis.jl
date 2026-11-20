# Heard et al. 2024: age structures and sample landscapes

- **Source**: Zenodo, "Data and code for 'Age truncation due to disease shrinks
  metapopulation viability for amphibians'", https://doi.org/10.5281/zenodo.10457229.
  Retrieved 2026-11-19. The full archive (R scripts, posterior samples, the compiled
  simulation outputs) is kept outside the package in `lab/papers/heard2024/zenodo/`.
- **Paper**: Heard GW, Scroggie MP, Hollanders M, Scheele BC (2024). Age truncation due to
  disease shrinks metapopulation viability for amphibians. *J Anim Ecol* 93:1670–1683.
  https://doi.org/10.1111/1365-2656.14177
- **Licence**: CC BY 4.0 (as labelled on the Zenodo record). Attribution: Heard et al. 2024,
  Zenodo.

## Files

| file | content |
|---|---|
| `age_structure.csv` | `Age_structure_data.csv`, unmodified: 403 frogs (*Litoria aurea* `Laur`, *L. raniformis* `Lran`), with `Epoch` (`Historical` = museum specimen collected before Bd emerged, `Contemporary`), `Site`, `Region`, coordinates, `SUL` (mm), `Sex`, `LAG` (lines of arrested growth = age in years) and `Reliability` of the reading. |
| `landscapes.csv` | the 20 bands of `sample_landscapes.grd/.gri` (25 × 25 cells of 65 m, UTM zone 55S, FLT4S, band-interleaved by line), converted to one row per landscape and grid row: `landscape,row,c1..c25`, 1 = habitat patch, 0 = not. Row 1 is the northern edge. |
| `landscape_connectivity.csv` | `sample_landscapes_connectivity.csv`, unmodified: connectivity index of each landscape. |
