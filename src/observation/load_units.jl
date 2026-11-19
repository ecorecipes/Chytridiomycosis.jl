# Load observation layer, part 1: the units map from a qPCR reaction to the load bases of
# src/units.jl.
#
#   Ct ──standard curve──▶ quantity per reaction ──extraction, dilution──▶ per swab
#      zoospore-count standards of strain s  ──▶ ZE (strain s)           ZEPerSwab
#      amplicon (ITS1 copy) standards         ──▶ ITS1 copies            ITS1CopiesPerSwab
#   ZE (strain s) ◀── × n_s ──▶ ITS1 copies ◀── × n_t ──▶ ZE (strain t)   restandardize
#
# n_s is the ITS1 copy number per zoospore of strain s (Longo et al. 2013, Fig. 4A, eight
# strains, 10 to 144). The sources give no conversion from ZE per swab to sporangia per host,
# to zoospores on the whole host, or from DNA to live zoospores beyond one lab study
# (Maguire et al. 2016); those steps are never taken silently (see swab_distribution).

const AVOGADRO = 6.02214076e23   # mol⁻¹, exact (SI 2019)
# Standard molecular weight of double-stranded DNA: (617.96 × bp + 36.04) g/mol, the formula
# of the common dsDNA mass calculators (anhydrous nucleotide pairs plus the 5′ and 3′ ends).
const DSDNA_G_PER_MOL_BP = 617.96
const DSDNA_G_PER_MOL_ENDS = 36.04

"""
    QPCRStandardCurve(slope, intercept, basis)
    QPCRStandardCurve(:longo2013, strain; standards = :zoospores)

A qPCR standard curve `Ct = slope · log10(q) + intercept`, where `q` is the quantity of
template in the reaction on the stated `basis`: `:zoospores` (zoospore genomic equivalents,
from counted zoospores) or `:dna_ng_per_uL` (equimolar genomic DNA).

The keyed form reads the curves of Longo et al. (2013, Table 1) for one of the eight
strains of [`longo_strains`](@ref): `standards = :zoospores` is Standard Set A and
`standards = :dna` is Standard Set B. See [`qpcr_efficiency`](@ref),
[`standard_quantity`](@ref) and [`cycle_threshold`](@ref).
"""
struct QPCRStandardCurve
    slope::Float64
    intercept::Float64
    basis::Symbol
    function QPCRStandardCurve(slope::Real, intercept::Real, basis::Symbol)
        slope < 0 || throw(ArgumentError("a standard curve has a negative slope, got $slope"))
        basis in (:zoospores, :dna_ng_per_uL, :copies) ||
            throw(ArgumentError("basis must be :zoospores, :dna_ng_per_uL or :copies"))
        return new(Float64(slope), Float64(intercept), basis)
    end
end

function QPCRStandardCurve(key::Symbol, strain::AbstractString; standards::Symbol = :zoospores)
    key === :longo2013 || throw(ArgumentError("only :longo2013 standard curves are in the compendium"))
    q = _longo_strain(strain)
    standards === :zoospores && return QPCRStandardCurve(q.ct_slope_A, q.ct_intercept_A, :zoospores)
    standards === :dna && return QPCRStandardCurve(q.ct_slope_B, q.ct_intercept_B, :dna_ng_per_uL)
    throw(ArgumentError("standards must be :zoospores (Set A) or :dna (Set B)"))
end

"""
    qpcr_efficiency(slope) -> Float64
    qpcr_efficiency(curve::QPCRStandardCurve)

Amplification efficiency `10^(−1/slope) − 1` of a standard curve with slope `slope` (Ct per
log10 of template): 1 means the product doubles every cycle (slope ≈ −3.32). Longo et al.
(2013, Table 1) print it as a percentage.
"""
qpcr_efficiency(slope::Real) = exp10(-1 / slope) - 1
qpcr_efficiency(c::QPCRStandardCurve) = qpcr_efficiency(c.slope)

"""
    standard_quantity(curve::QPCRStandardCurve, ct) -> Float64

The template quantity per reaction that gives cycle threshold `ct` on `curve`, in the
curve's basis: `10^((ct − intercept)/slope)`.
"""
standard_quantity(c::QPCRStandardCurve, ct::Real) = exp10((ct - c.intercept) / c.slope)

"""
    cycle_threshold(curve::QPCRStandardCurve, q) -> Float64

The cycle threshold of `q` units of template on `curve`; the inverse of
[`standard_quantity`](@ref).
"""
function cycle_threshold(c::QPCRStandardCurve, q::Real)
    q > 0 || throw(DomainError(q, "a template quantity must be positive"))
    return c.slope * log10(q) + c.intercept
end

"""
    longo_strains() -> Vector{String}

The eight Bd strains of Longo et al. (2013) in the order of their Table 1: JEL404, MexMkt,
JEL427, PAB01, LFT001_01, CLFT023, CLFT024, LBAbercrom. Each is a named set of
`data/compendium/longo2013.toml`.
"""
function longo_strains()
    sets = String[]
    for r in parameter_records(:longo2013)
        isempty(r.set) || r.set in sets || push!(sets, r.set)
    end
    return sets
end

function _longo_strain(strain::AbstractString)
    strain in longo_strains() ||
        throw(ArgumentError("unknown strain \"$strain\"; Longo et al. 2013 have $(longo_strains())"))
    return parameter_set(:longo2013; set = strain)
end

"""
    strain_copy_number(strain) -> Float64

ITS1 (ITS1/5.8S) copies per zoospore of a Longo et al. (2013) strain, measured against an
amplicon standard (Fig. 4A): from 10 (MexMkt) to 144 (CLFT023). This is the
`copies_per_zoospore` of [`convert_load`](@ref) for loads quantified with that strain's
zoospore standards.
"""
strain_copy_number(strain::AbstractString) = _longo_strain(strain).copies_per_zoospore

_copy_number(n::Real) = n > 0 ? Float64(n) :
    throw(ArgumentError("a copy number per zoospore must be positive, got $n"))
_copy_number(s::AbstractString) = strain_copy_number(s)

"""
    restandardize(x::ZEPerSwab, from => to) -> ZEPerSwab

Re-express a load in zoospore equivalents (ZE) measured against the zoospore standards of
one strain as ZE of another. `from` and `to` are Longo et al. (2013) strain names, or ITS1
copy numbers per zoospore. The swab holds `x × n_from` ITS1 copies, which is
`x × n_from / n_to` ZE of the other strain: two explicit [`convert_load`](@ref) steps through
[`ITS1CopiesPerSwab`](@ref).

For example, a swab read as 1,000 ZE against JEL404 standards (39 copies per zoospore) is
271 ZE against CLFT023 standards (144 copies), and 3,900 ZE against MexMkt standards (10).
It assumes, like the standards themselves, that the swab's strain has the standard's copy
number; when the infecting strain differs, ZE count copies, not zoospores.
"""
function restandardize(x::ZEPerSwab, strains::Pair)
    nf, nt = _copy_number(strains.first), _copy_number(strains.second)
    copies = convert_load(x, ITS1CopiesPerSwab; copies_per_zoospore = nf)
    return convert_load(copies, ZEPerSwab; copies_per_zoospore = nt)
end

"""
    dsdna_mass(bp; bp_mass = nothing) -> Float64

Mass in ng of one molecule of double-stranded DNA of `bp` base pairs. By default the
standard formula `(617.96 bp + 36.04)` g/mol; with `bp_mass` (Da per base pair) it is
`bp × bp_mass` g/mol, the form Longo et al. (2013) print (`bp_mass = 660`). Divided by
Avogadro's number and scaled to ng.

The copy number of an amplicon standard is its DNA mass divided by this, so the choice of
formula scales every copy number derived from it.
"""
function dsdna_mass(bp::Real; bp_mass::Union{Nothing,Real} = nothing)
    bp > 0 || throw(ArgumentError("bp must be positive"))
    gmol = bp_mass === nothing ? DSDNA_G_PER_MOL_BP * bp + DSDNA_G_PER_MOL_ENDS : bp * bp_mass
    return gmol / AVOGADRO * 1e9
end

"""
    genome_equivalents(ng; genome_bp, bp_mass = nothing) -> Float64

The number of haploid-set genomes (zoospores, for Bd) whose DNA weighs `ng` nanograms, for a
genome of `genome_bp` base pairs per zoospore ([`dsdna_mass`](@ref) gives the mass of one).
Longo et al. (2013) use a 46 Mb diploid genome and 660 Da per base pair:
`genome_equivalents(0.3; genome_bp = 46e6, bp_mass = 660)` is about 5,950, their "~6,000
zoospores" per µL of the top DNA standard.
"""
genome_equivalents(ng::Real; genome_bp::Real, bp_mass::Union{Nothing,Real} = nothing) =
    ng / dsdna_mass(genome_bp; bp_mass)

"""
    shed_viability(day = 0; key = :maguire2016, reference = true) -> Float64

The fraction of zoospores shed from a live frog that are alive (stain excludes propidium
iodide), from the mixed logistic regression of Maguire et al. (2016, Table 2):
`logit v = b0 + b_day · day` for the reference species, plus the species effect and the
day × species interaction with `reference = false`. Coefficients from the compendium.

qPCR counts the DNA of dead zoospores too. A model whose pool holds live zoospores predicts a
qPCR (DNA) count of `live / v`; [`qpcr_zoospores`](@ref) and [`viable_fraction`](@ref) do this
for the zoospore-pool models.

The table does not name the reference species. The text says the odds of viability rose for
*Atelopus zeteki* and not for *Hylomantis lemur*, which fits *A. zeteki* as the reference. On
day 0 it is 0.78, the "approximately 80%" of the abstract. One lab study of one strain (JEL
423): not a general constant.
"""
function shed_viability(day::Real = 0; key::Symbol = :maguire2016, reference::Bool = true)
    q = parameter_set(key)
    z = q.viability_logit_intercept + q.viability_logit_day * day
    reference || (z += q.viability_logit_species + q.viability_logit_day_species * day)
    return 1 / (1 + exp(-z))
end

# Validation helpers ------------------------------------------------------------------------

# Number of Longo strains whose printed efficiency (Table 1, Standard Set :A or :B) follows
# from the printed slope as 10^(−1/slope) − 1, allowing for the rounding of both: the slope
# to 0.01 and the efficiency to its printed digits (0.5 for whole percentages, else 0.05).
function _longo_efficiency_consistent(set::Symbol)
    n = 0
    for s in longo_strains()
        q = parameter_set(:longo2013; set = s)
        slope, eff = set === :A ? (q.ct_slope_A, q.efficiency_A) : (q.ct_slope_B, q.efficiency_B)
        lo, hi = extrema(100 .* qpcr_efficiency.((slope - 0.005, slope + 0.005)))
        r = isinteger(eff) ? 0.5 : 0.05
        n += (lo - r <= eff <= hi + r)
    end
    return n
end

# Whether one strain's zoospore-count standard curve (Set A) has the highest Ct of all eight
# at every quantity on a grid over the standards' range (0.1 to 1000 zoospores).
function _longo_highest_ct(strain::AbstractString)
    curves = Dict(s => QPCRStandardCurve(:longo2013, s) for s in longo_strains())
    return all(exp10.(range(-1, 3; length = 41))) do q
        c = cycle_threshold(curves[strain], q)
        all(cycle_threshold(curves[s], q) < c for s in longo_strains() if s != strain)
    end
end
