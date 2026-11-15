# Hollanders et al. 2023: *Mixophyes fleayi* robust-design capture–recapture data

- **Paper**: Hollanders M, Grogan LF, Nock CJ, McCallum HI, Newell DA (2023). Recovered frog
  populations coexist with endemic *Batrachochytrium dendrobatidis* despite load-dependent
  mortality. *Ecological Applications* 33:e2724. https://doi.org/10.1002/eap.2724 (open
  access, CC BY 4.0).
- **Data**: Dryad, https://doi.org/10.5061/dryad.g1jwstqtb, version 4 of 2022-08-03.
  **Licence: CC0 1.0**, from the Dryad API record on 2026-11-15.
- **How the files were obtained.** On 2026-11-15 the Dryad file endpoints refused the
  automated client (HTTP 401/403). The same files were downloaded from the authors' public
  repository `mhollanders/mfleayi-adults` at commit
  `a325a4c61a5687c25d7ccc6b30d23a1c25ecc390`, under `data/`. Every file here matches,
  byte for byte, the SHA-256 digest that the Dryad API lists for its Dryad copy:

  | file | SHA-256 |
  |---|---|
  | `bw.csv` | `2e382b08e80c57d7608999c469e92fc360086f196f0e22598330a2168c623ee8` |
  | `ch.csv` | `08684be90fed9dfae27fdc55b01629df0b1bd123c26c2493e8f0422ccebce334` |
  | `covs-int.csv` | `e62335cf47288d6befe7ecb6c61b96af66bef86d6ec47dbf423a4ca52e03ffbf` |
  | `covs-surv.csv` | `121f601ccef0240d335475bc4c69fc5dc6dcfd20d26b63ba260ec822f0f6e198` |
  | `dates.csv` | `03fb0cf6908ba7067445e19f48b8b02c2d379942038ef467055168704190584d` |
  | `ints.csv` | `12d339b48f490b3958059ab8f69d7ad163bb09abd1ddf570e18a359afefa2fd2` |
  | `length.csv` | `3814701a9a4ec4b8487cd502b2ec52f4f48cdbf17e44b8ea4cebb2658b5414e1` |
  | `loads.csv` | `222db1d45082861d7fb0ee05a7a1c1f124ce061801522c28c533334ed77b69a1` |
  | `sex.csv` | `8624499e2a64acb72748b420d66d8cc916388b35f2f26e67c739eb1c9907d44c` |

  `community.csv`, the swabs of sympatric species for the paper's first model, was not
  bundled: that model is not implemented. The full download, including that file, is in
  `lab/papers/hollanders2023_si/`.
- **Code**: Zenodo, https://doi.org/10.5281/zenodo.6981761 (MIT licence; `m1-infections.qmd`
  and `m2-cmr.qmd`, retrieved 2026-11-15). It was read, not run.

## Layout (from `m2-cmr.qmd`)

The CSVs hold R arrays written column by column. [`load_hollanders2023`](../../src/models/hollanders2023.jl)
reshapes them the same way.

| file | shape | content |
|---|---|---|
| `ch.csv` | 462 × (23 · 3 · 3) | capture code by individual, then primary (fastest), secondary and site: 0 not seen, 1 seen *Bd*−, 2 seen *Bd*+, NA no survey or padding |
| `loads.csv` | 462 × (23 · 3) | swab load, ITS1 copies per swab, by primary within site; NA when not captured or not infected |
| `bw.csv`, `length.csv` | 462 × (23 · 3) | body mass (g) and snout–urostyle length (mm) |
| `sex.csv` | 462 × 3 | 1 female, 0 male, NA unknown or padding |
| `ints.csv` | 22 × 3 | interval lengths between primaries, in six-week units |
| `dates.csv` | 23 × 3 | primary dates (d/m/yy); Tuntable primary 22 reads `1/2/20` for 1 February 2021 |
| `covs-int.csv` | 69 rows | SILO daily maximum and minimum temperature (°C) and rainfall (mm) over the interval ending at each primary, site by site |
| `covs-surv.csv` | 207 rows | SILO maximum and minimum temperature, vapour pressure and daily rainfall on each secondary survey day |

Sites, in order: Brindle Creek (462 individuals), Tuntable Creek (136) and Bat Cave Creek (88).

## The authors' posterior draws

`mcmc-samples/m2-draws.csv` (54 MB, 10 000 draws) comes from the same repository. The
repository has no licence, so the draws are **not bundled**.
[`hollanders2023_draws_path()`](../../src/models/hollanders2023.jl) finds them at
`ENV["CHYTRIDIOMYCOSIS_HOLLANDERS2023_DRAWS"]`, or at
`papers/hollanders2023_si/repo_mfleayi-adults/mcmc-samples/m2-draws.csv` under any ancestor
of the package directory. The tests and validation targets that use the draws are skipped
when the file is absent.
