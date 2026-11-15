# Palomar et al. 2023: *Bufo spinosus* multi-event capture–recapture data (not available)

- **Paper**: Palomar G, Fernández-Chacón A, Bosch J (2023). Amphibian survival compromised by
  long-term effects of chytrid fungus. *Biodiversity and Conservation* 32:793–809.
  https://doi.org/10.1007/s10531-022-02525-3 (open access, CC BY 4.0).
- **Data statement**: "Data will be available online upon publication."
- **Search, 2026-11-15**: no deposit was found. The following were searched:
  - Dryad, Zenodo and DataCite, for "Bufo spinosus" with survival or chytrid, and for
    Palomar with Bosch;
  - figshare;
  - the DIGITAL.CSIC repository record found by search, which turned out to be a different
    paper (Thumsová et al. 2023, *Oecologia*).

  The supplementary material (Online Appendices A and B) gives only the model pattern and
  the goodness-of-fit table. It is in `lab/papers/palomar2023_si1.docx` and `_si2.docx`.
- **Consequence**: nothing is bundled and nothing is refitted.
  [`palomar2023_structure`](../../src/models/palomar2023.jl) and
  [`palomar2023_model_selection`](../../src/fitting/palomar2023.jl) take histories in the
  appendix coding:
  - 0: not seen;
  - 1: seen, infected;
  - 2: seen, uninfected;
  - 3: seen, status unknown.

  They will fit the 13 models of Table 1 once the 632 male histories are available. The
  published estimates are registered as `pending` validation targets.
