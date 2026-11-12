# Units layer (lab/wiki/parameters/units.md): typed loads, explicit conversions, time units.
@testset "Units (Phase 0)" begin
    # Loads store log10 and read back exactly.
    z = ZEPerSwab(1e4)
    @test log10(z) == 4 && load_value(z) ≈ 1e4
    @test ZEPerSwab(; log10 = 4) == z
    @test load_value(ZEPerSwab(0)) == 0
    @test_throws DomainError ZEPerSwab(-1)
    @test sprint(show, z) == "ZEPerSwab(10000.0)"
    @test ZEPerSwab(10) < ZEPerSwab(100) && ZEPerSwab(100) ≈ ZEPerSwab(100.0000001)

    # Round trips, with the assumptions stated.
    q = parameter_set(:longo2013)
    @test (q.copies_per_zoospore_min, q.copies_per_zoospore_max) == (10.0, 144.0)
    c = ITS1CopiesPerSwab(5000)
    ze = convert_load(c, ZEPerSwab; copies_per_zoospore = q.copies_per_zoospore_min)
    @test load_value(ze) ≈ 500
    @test convert_load(ze, ITS1CopiesPerSwab; copies_per_zoospore = q.copies_per_zoospore_min) ≈ c
    zx = convert_load(c, ZEPerSwab; copies_per_zoospore = 50, extraction_factor = 10)
    @test load_value(zx) ≈ 1000
    @test convert_load(zx, ITS1CopiesPerSwab; copies_per_zoospore = 50, extraction_factor = 10) ≈ c
    # The strain's copy number alone spans more than a log: 10 to 144 copies per zoospore.
    lo = convert_load(c, ZEPerSwab; copies_per_zoospore = q.copies_per_zoospore_max)
    @test log10(ze) - log10(lo) ≈ log10(14.4)
    on_host = convert_load(z, Zoospores; swab_fraction = 0.1)
    @test load_value(on_host) ≈ 1e5 && convert_load(on_host, ZEPerSwab; swab_fraction = 0.1) ≈ z
    @test convert_load(z, ZEPerSwab) === z
    # Assumptions are required, not defaulted.
    @test_throws UndefKeywordError convert_load(c, ZEPerSwab)
    @test_throws ArgumentError convert_load(c, ZEPerSwab; copies_per_zoospore = 0)
    @test_throws ArgumentError convert_load(z, Zoospores; swab_fraction = 1.5)
    @test_throws ArgumentError convert_load(c, SporangiaPerHost; ze_per_sporangium = 1.0)

    # ZE to sporangia has no published basis: it warns.
    s = @test_logs (:warn, r"no published conversion") convert_load(z, SporangiaPerHost;
        ze_per_sporangium = 2.0)
    @test load_value(s) ≈ 5e3
    @test (@test_logs (:warn,) convert_load(s, ZEPerSwab; ze_per_sporangium = 2.0)) ≈ z

    # Briggs Smax (sporangia per host) and the 10,000 ZE crash load are different types.
    smax = parameter_load(:briggs2010, :Smax)
    crash = parameter_load(:kinney2011, :crash_load)
    @test smax isa SporangiaPerHost && crash isa ZEPerSwab
    @test load_value(smax) == load_value(crash) == 1e4
    @test_throws ArgumentError smax == crash
    @test_throws ArgumentError smax < crash
    @test_throws ArgumentError isapprox(smax, crash)
    @test_throws ArgumentError parameter_load(:briggs2010, :γ)          # not a load
    @test_throws ArgumentError parameter_load(:briggs2010, :no_such_symbol)

    # Rates and durations in explicit time units.
    @test convert_rate(52.0, :year => :day; days_per_year = 364) ≈ 1 / 7
    @test convert_rate(1.0, :day => :hour) ≈ 1 / 24
    @test convert_rate(convert_rate(0.3, :week => :year), :year => :week) ≈ 0.3
    @test_throws ArgumentError convert_rate(1.0, :day => :fortnight)
    recs = parameter_records(:schmidt2017)
    e = recs[findfirst(r -> r.symbol === :e, recs)]                    # 364 yr⁻¹
    @test in_time_unit(e, :day; days_per_year = 364) ≈ 1.0
    β = recs[findfirst(r -> r.symbol === :β, recs)]                    # "ha yr^-1"
    @test in_time_unit(β, :week; days_per_year = 364) ≈ 1.0
    @test in_time_unit(e, :year) == e.value
    wd = parameter_records(:woodhams2008)
    s23 = wd[findfirst(r -> r.symbol === :s && r.set == "23C", wd)]    # duration, hours
    @test in_time_unit(s23, :day) ≈ 87.7 / 24
    K = recs[findfirst(r -> r.symbol === :K, recs)]
    @test_throws ArgumentError in_time_unit(K, :day)
    @test survival_to_hazard(hazard_to_survival(0.2, 3.0), 3.0) ≈ 0.2
    @test survival_to_hazard(0.5, 7.0) ≈ log(2) / 7
    @test_throws DomainError survival_to_hazard(0.0, 1.0)

    # Declared model time units.
    @test time_unit(SchmidtParameters(:schmidt2017)) === :year
    @test time_unit(BriggsParameters(:briggs2010)) === :day
    @test time_unit(WoodhamsParameters) === :hour
    @test_throws ArgumentError time_unit(1.0)
end
