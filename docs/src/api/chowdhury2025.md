# Life stages and temperature (Chowdhury 2025)

```@meta
CurrentModule = Chytridiomycosis
```

Chowdhury et al. (2025, *Ecol Model* 501:111028,
[doi:10.1016/j.ecolmodel.2025.111028](https://doi.org/10.1016/j.ecolmodel.2025.111028))
combined laboratory experiments on the eastern newt (*Notophthalmus viridescens*) with a
compartmental model of *B. salamandrivorans* (Bsal). They asked how temperature (6 and
14 °C) and host density change epidemics in a host whose life stages live in different
habitats. Adults and larvae live in water and efts (juveniles) on land. The model was
checked against the published PDF and its online supplement (Appendices A–H). No code or
data were deposited.

## The model

Larvae `S_l` are never infected. They metamorphose at rate `η_j`, and a fraction `α`
survives to become efts. Efts and adults each pass through susceptible, exposed, infectious
and recovered classes (`S, E, I, R`). Efts join the adults at rate `η_a(1 − A/K)`, which
caps the adult density near `K`. The rest of the efts that mature leave the system. Each
habitat has its own zoospore pool: `Z_j` in soil and `Z_a` in water. Time is in days and
densities are per m²:

```math
\begin{aligned}
\dot S_l &= \lambda A - \eta_j S_l,\\
\dot S_j &= \alpha\eta_j S_l - \beta_j \frac{S_j I_j}{N_T} - f_j(Z_j) S_j - (\delta_j + \eta_a) S_j,\\
\dot E_j &= f_j(Z_j) S_j + \beta_j \frac{S_j I_j}{N_T} - (\sigma_j + \delta_j + \eta_a) E_j,\\
\dot I_j &= \sigma_j E_j - (\gamma_j + \delta_j + \mu_j + \eta_a) I_j,\\
\dot R_j &= \gamma_j I_j - (\delta_j + \eta_a) R_j,\\
\dot S_a &= \eta_a\Big(1 - \frac{A}{K}\Big) S_j - \beta_a(A) \frac{S_a I_a}{N_A} - f_a(Z_a) S_a - \delta_a S_a,\\
\dot E_a &= f_a(Z_a) S_a + \eta_a\Big(1 - \frac{A}{K}\Big) E_j + \beta_a(A) \frac{S_a I_a}{N_A} - (\sigma_a + \delta_a) E_a,\\
\dot I_a &= \eta_a\Big(1 - \frac{A}{K}\Big) I_j + \sigma_a E_a - (\gamma_a + \delta_a + \mu_a) I_a,\\
\dot R_a &= \eta_a\Big(1 - \frac{A}{K}\Big) R_j + \gamma_a I_a - \delta_a R_a,\\
\dot Z_j &= \omega_j I_j - \xi_j Z_j, \qquad \dot Z_a = \omega_a I_a - \xi_a Z_a.
\end{aligned}
```

Here `A` is the number of adults, `N_T` the number of efts and `N_A = A + S_l` the aquatic
population. Infection comes by two routes:

- **Direct contact.** Adult contacts saturate with density, `β_a(A) = P_a c_a A/(χ_a + A)`
  (Holling type II, eq. 3a). Eft contacts are rare and constant, `β_j = c_j P_j` (eq. 3b).
- **Zoospores.** `f_i(Z) = ρ_i b_i Z/(κ_i + Z)` (eq. 4). `ρ` is the fraction of contacted
  spores that encyst, `b` a contact rate with spores, and `κ` the ID50.

[`BsalStageTemperatureProblem`](@ref) lowers the model to an `ODEProblem`. The direct terms
are [transmission functions](bsal_transmission.md). The adult default is
`HollingContact(P_a, c_a, χ_a)` diluted by larvae (`A/N_A`), and the eft default is
`FrequencyDependent(c_j P_j)`. Swap either to test another contact structure.

## What the paper does not say, and where it disagrees with itself

- **`χ_a` is not printed.** At 14 °C the Fig. 1 curve, the starting `β_a` of Fig. 7 and the
  `R_a` column of Table 3 all give `χ_a ≈ 5.9` per m². At 6 °C, Fig. 7 with Table 2's
  `c_a = 464.77` gives `χ_a ≈ 12`, and these values reproduce Figs 4, 5 and 7. But the
  Fig. 1 curve at 6 °C and the 6 °C `R_a` column of Table 3 need `c_a ≈ 729` and
  `χ_a ≈ 23.3`. So Table 2's `c_a` at 6 °C was not used for Table 3.
- **Table 3 uses a different eft contact rate.** §2.2 counts 8 contacts in 14 days among 12
  efts and sets `c_j = 8/14 = 0.571`. That is per enclosure, not per eft. Table 3's `R_j`
  needs `8/(14 × 12) = 0.0476`, and also `b_j = 0.005` (§6.2). With these, 14 °C
  reproduces, but at 6 °C the environmental part of `R_j` is three times the formula's.
- **Spore persistence in soil.** §2.3 says 25 days at 6 °C and 20 at 14 °C. Table 2 gives
  the reverse (`ξ_j = 1/20, 1/25`), and only Table 2 reproduces Table 3. Appendix C says
  efts mature after 3.5 years and Table 2 says 3; again only Table 2 reproduces Table 3.
  `as_printed = true` uses the text values.
- **Efts.** With Table 2's birth rate (`λ = 0.82` eggs per adult per day) larvae recruit
  several efts per m² within months, so the simulated eft numbers rise. In Figs 4, 5 and
  E.1 they stay flat or fall and never recover over 10 years. With negligible recruitment
  (`λ ≈ 0`) the 14 °C eft curves are close to the paper's. The adult curves are the same
  either way.
- **Units of the pools.** `Z_a` is compared with `κ_a` in spores per mL and `Z_j` with `κ_j`
  in spores per m², but both pools are filled by `ω I` with `I` per m². The paper equates
  them without a volume. The package keeps the paper's numbers.

`ChowdhuryParameters(:chowdhury2025; variant = :table3)` gives the values that reproduce
Table 3.

## What reproduces

| Result | Paper | Here |
|---|---|---|
| Adult infected peak, 16 m⁻² (Fig. 4) | ≈ 8.3 (6 °C), 7.5 (14 °C) | 8.3, 7.5 |
| Aquatic zoospore peak, 16 m⁻² (Fig. 5) | ≈ 8,300 (6 °C), 17,900 (14 °C) | 8,300, 17,900 |
| ID50 exceeded by `Z_a` (Fig. 5) | 16 m⁻² only at 6 °C; 3 of 4 at 14 °C | same |
| `β_a` at the start (Fig. 7) | ≈ 51 … 12 (6 °C), 218 … 75 (14 °C) | same (by construction of `χ_a`) |
| `R_a`, Table 3, 14 °C | 16.55, 27.85, 43.72, 65.67 | 16.55, 27.86, 43.75, 65.73 |
| `R_j`, Table 3, 14 °C (`:table3`) | 1.26, 2.22, 4.15, 8.01 | 1.26, 2.22, 4.14, 7.99 |
| `R_a`, Table 3, 6 °C (`:table3`) | 4.82, 9.1, 16.5, 28.1 | 4.83, 9.11, 16.48, 28.11 |
| `R_j`, Table 3, 6 °C | 2.35, 4.33, 8.29, 16.21 | not reproduced (1.03 … 5.66) |

The validation targets (`validation_targets(; source = :chowdhury2025)`) record each claim.
Table 3's numbers mix parameter sets: no single parameter set reproduces Figs 4–7 and
Table 3 at 6 °C.

## Reproduction numbers

The next-generation matrix has two non-zero eigenvalues (Appendix F), one per habitat. Each
is the sum of a direct and an environmental route:

```math
\mathcal{R}_j = \frac{c_jP_j\sigma_j}{D_j(\eta_a + \delta_j + \sigma_j)}
 + \frac{\rho_j S_j^*\omega_j\sigma_j b_j}{\kappa_j\xi_j D_j(\eta_a + \delta_j + \sigma_j)},
\qquad
\mathcal{R}_a = \frac{S_a^{*2}c_aP_a\sigma_a}{(\sigma_a + \delta_a)D_a(S_l^* + S_a^*)(\chi_a + S_a^*)}
 + \frac{\rho_a\omega_a b_a\sigma_a S_a^*}{\kappa_a\xi_a(\sigma_a + \delta_a)D_a},
```

with `D_j = δ_j + η_a + γ_j + μ_j` and `D_a = μ_a + δ_a + γ_a`. `R₀ = max(R_j, R_a)`. The
package computes them from the numerical next-generation matrix, so they also hold for
other contact functions, and checks them against these formulas in the tests. The
disease-free state is far from the simulated one: it has about 16 efts and 85 larvae per
adult.

## Sensitivity

The paper ranks parameters by Latin hypercube sampling with partial rank correlation
(Fig. 8, Appendix H) and by elasticities of `R_j` and `R_a` (Appendix G). The sampling
ranges are not given, so the PRCC values cannot be reproduced. [`stage_elasticities`](@ref)
computes the elasticities, and their signs agree with every direction stated in §6.1.

## API

```@docs
ChowdhuryParameters
BsalStageTemperatureProblem
chowdhury_scenarios
stage_transmission_rates
stage_totals
stage_decline_time
stage_disease_free_state
stage_next_generation_matrix
stage_reproduction_numbers
stage_elasticities
```

`basic_reproduction_number(p::ChowdhuryParameters)` returns `max(R_j, R_a)`.
