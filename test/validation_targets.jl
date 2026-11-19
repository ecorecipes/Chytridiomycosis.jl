# Validation-target registry (data/validation/targets.toml).
# Runs every registered target that has an implementation. Discrepant targets are
# @test_broken against the paper and @test against our recorded value.
@testset "Validation targets (Phase 0)" begin
    ts = validation_targets()
    @test allunique(t.id for t in ts)
    @test all(t -> t.status in (:reproduced, :qualitative, :discrepant, :pending), ts)
    @test validation_target("schmidt2017_K_threshold").source === :schmidt2017
    @test_throws ArgumentError validation_target("no_such_target")
    @test_throws ArgumentError reproduce(validation_target("blooi2015_heat_clearance"))
    @test length(validation_targets(; status = :discrepant)) == 63
    @test all(t -> t.source === :schmidt2017, validation_targets(; source = :schmidt2017))
    # Wilber 2022 targets need the authors' deposit, which is not bundled (data/wilber2022)
    w22 = Chytridiomycosis._wilber2022_find() !== nothing
    # Hollanders 2023 targets that read the authors' posterior draws (not bundled)
    h23 = Chytridiomycosis._hollanders2023_draws_find() !== nothing
    h23draws(t) = occursin(r"hollanders2023_posterior|_hollanders2023_slope_deviation", t.call)
    @testset "$(t.id)" for t in filter(t -> t.status !== :pending &&
                                             (w22 || t.source !== :wilber2022) &&
                                             (h23 || !h23draws(t)), ts)
        r = check_target(t)
        if t.status === :discrepant
            @test_broken r.published
            @test r.regression
        else
            @test r.published
            r.regression === nothing || @test r.regression
        end
    end
end
