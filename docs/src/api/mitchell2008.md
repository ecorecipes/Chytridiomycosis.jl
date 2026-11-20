# Bd in common toads at a pond (Mitchell et al. 2008)

```@meta
CurrentModule = Chytridiomycosis
```

This is a reference implementation of the pond model of Mitchell et al. (2008), *Proc R Soc
B* 275:329–334 ([doi:10.1098/rspb.2007.1356](https://doi.org/10.1098/rspb.2007.1356)): common
toads (*Bufo bufo*) breeding at a pond, with life stages and Bd zoospores that persist in the
water outside the host. It is checked against the paper. Parameters are read from the
[compendium](compendium.md): `MitchellParameters(:mitchell2008)`.

```@docs
MitchellParameters
MitchellPondProblem
```
