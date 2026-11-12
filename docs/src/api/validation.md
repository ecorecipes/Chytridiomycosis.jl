# Validation targets

```@meta
CurrentModule = Chytridiomycosis
```

The registry in `data/validation/targets.toml` lists the published numbers and outcomes that
the implementations must reproduce. Each target has:

- the paper and the place in it;
- the package call that reproduces it;
- a tolerance;
- a status: `reproduced`, `qualitative`, `discrepant` (the paper's number is not reproduced,
  and our value is recorded as a regression check) or `pending` (no implementation yet).

The test suite runs every target that is not pending.

```julia
for t in validation_targets(; status = :discrepant)
    println(t.id, ": paper ", t.published, ", we get ", check_target(t).value)
end
```

```@docs
ValidationTarget
validation_targets
validation_target
reproduce
check_target
```
