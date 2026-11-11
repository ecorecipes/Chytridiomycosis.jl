# Parameter compendium

```@meta
CurrentModule = Chytridiomycosis
```

Published parameter values are read from TOML files in `data/compendium/`, one per source
paper. Each record has its units, a provenance code (`M` measured, `F` fitted, `E` elicited,
`A` assumed, `D` derived), a role, the system it applies to, where in the paper it appears,
and whether it has been checked against the paper (`verified`). Models take their published
values only from here.

```@docs
compendium_sources
parameter_records
parameter_set
ParameterRecord
```
