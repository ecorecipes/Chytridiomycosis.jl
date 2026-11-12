# Validation-target registry (data/validation/targets.toml; lab/wiki/validation/targets.md).
# Runs every registered target that has an implementation. Discrepant targets are
# @test_broken against the paper and @test against our recorded value.
@testset "Validation targets (Phase 0)" begin
    ts = validation_targets()
    @test allunique(t.id for t in ts)
    @test all(t -> t.status in (:reproduced, :qualitative, :discrepant, :pending), ts)
    @test validation_target("schmidt2017_K_threshold").source === :schmidt2017
    @test_throws ArgumentError validation_target("no_such_target")
    @test_throws ArgumentError reproduce(validation_target("blooi2015_heat_clearance"))
    @test length(validation_targets(; status = :discrepant)) == 5
    @test all(t -> t.source === :schmidt2017, validation_targets(; source = :schmidt2017))
    @testset "$(t.id)" for t in filter(t -> t.status !== :pending, ts)
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
