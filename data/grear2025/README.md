# Grear et al. 2025: US pan-amphibian Bd swab survey with temperature covariates

- **Paper**: Grear DA, Adams MJ, Backlin AR, et al. (2025). Pan-amphibia distribution of the
  fungal parasite *Batrachochytrium dendrobatidis* varies with species and temperature.
  *Ecological Monographs* 95:e70001. https://doi.org/10.1002/ecm.70001 (open access,
  CC BY-NC 4.0). The full text could not be fetched by an automated client (the publisher's
  bot check); only the abstract (USGS Publications Warehouse, record 70263881) was read.
- **Source of this file**: the authors' USGS software release, Grear DA (2025)
  "Pan-Amphibia Bd distribution: software release", https://doi.org/10.5066/P1MXUJZM
  (repository https://code.usgs.gov/nwhc-software/pan-amphibia-bd-distribution, commit
  `b279e8f08cbde26d03ddb63b35293cd2b12edbcd`, cloned 2026-11-19). Its README describes
  `BsalBdDataFull.csv` as "the dataset that will be posted to the associated data release";
  the code downloads the same file from the ScienceBase data release of Waddle et al. (2019),
  "Data from a national survey for the amphibian chytrid fungus *Batrachochytrium
  salamandrivorans*", https://doi.org/10.5066/P9BGQA1T, which has the covariates joined.
- **Licence**: public domain. The software release is marked CC0 1.0 ("Unless otherwise
  noted, this project is in the public domain in the United States because it contains
  materials that originally came from the United States Geological Survey"); the ScienceBase
  release is a USGS data release (public domain). The temperature and precipitation
  covariates were computed by the authors from PRISM climate data.
- **What is bundled**: `grear2025_swabs.csv` keeps 15 of the 71 columns of
  `BsalBdDataFull.csv` (SHA-256
  `e965ccbc8e914870fce53857174fdbb9ca9d2010a8350f0af0f505598bbe198c`), all 11,189 rows, values
  unchanged (`NA` marks missing):

| column | meaning |
|---|---|
| `SwabID` | swab identifier |
| `CollectionDate` | m/d/yyyy |
| `Latitude`, `Longitude`, `State` | location |
| `Gesp`, `order` | species (genus species; `Taricha sp.` unresolved) and order (Caudata, Anura) |
| `Age`, `Sex` | Adult/Juvenile/Larval; Female/Male/Unknown |
| `siteID`, `reg` | site created by the authors from pairwise distances; NOAA climate region |
| `Bd.Result` | qPCR result as recorded |
| `copyNum` | Bd quantity per swab, the response of the authors' hurdle model (`0` = not detected) |
| `m30max` | mean daily maximum air temperature over the 30 days before sampling (°C) |
| `t10precip` | total precipitation over the 10 days before sampling (mm) |

Not bundled: the phylogeny `outputSppInData.nex` of the software release (one tree pruned
from the Jetz & Pyron 2018 amphibian posterior on vertlife.org, whose terms are not stated in
the release). `grear2025_phylogeny()` reads it from `ENV["CHYTRIDIOMYCOSIS_GREAR2025"]` or
`lab/papers/grear2025/code/src/` when present.
