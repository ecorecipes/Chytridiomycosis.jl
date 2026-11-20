# Sheets et al. 2021 Front Vet Sci 8:687084 (doi:10.3389/fvets.2021.687084): logistic growth of
# normalised MTT viability for five Bd isolates at seven temperatures. Oracles:
# - Table 2 (r and K of each isolate and temperature; mean viability) and Table 4;
# - the authors' R Markdown (figshare doi:10.6084/m9.figshare.14714865): its likelihood,
#   written out again here, and a one-off run of its fitlogistic() loop (set.seed(2021)) on
#   the deposited data, kept in lab/papers/sheets2021/oracle/logisticfits_R.csv; the values
#   used below are copied from that run;
# - the Fig. 4 legend (VT 26 °C: r = 8.58), which the R run reproduces;
# - the Supplementary Material (Data_Sheet_1.PDF, the rendered R Markdown): its printed fits
#   (compendium r_fit, K_fit) and profile-likelihood 95% intervals. It has no Johnson-Lewin fits.
using Chytridiomycosis
using Test
using LogDensityProblems
using Distributions: Normal, logpdf

@testset "Sheets 2021 Bd isolate thermal performance" begin
    C = Chytridiomycosis
    data = load_sheets2021()
    table(iso, T) = parameter_set(:sheets2021; set = "$(iso)_$(Int(T))C")

    @testset "data and likelihood" begin
        @test length(data) == 5509
        @test Set(d.isolate for d in data) == Set(C.SHEETS_ISOLATES)
        @test sort(unique(d.T for d in data)) == collect(C.SHEETS_TEMPERATURES)
        m = SheetsGrowthModel("LA", 21; data)
        @test C.nobs(m) == 38                                    # R: nobs = 38
        @test LogDensityProblems.dimension(m) == 5
        p = (x0 = 0.005, r = 1.5, K = 1.0, b0 = 0.01, b1 = 0.2)
        x = log.(collect(p))
        # obj() of the R code
        ll = sum(logpdf(Normal(μ, p.b0 + 1e-4 + p.b1 * μ), y)
                 for (t, y) in zip(m.t, m.y)
                 for μ in (p.K * p.x0 / (p.x0 + (p.K - p.x0) * exp(-p.r * t)),))
        @test LogDensityProblems.logdensity(m, x) ≈ ll
        @test_throws ArgumentError SheetsGrowthModel("CA", 21; data)
    end

    @testset "Table 2: logistic fits" begin
        # reproduced to the printed 4 significant figures (22 of the 35 fits)
        same = [("LA", 4), ("LA", 12), ("LA", 21), ("LA", 25), ("TN", 12), ("TN", 17),
                ("TN", 21), ("TN", 25), ("TN", 27), ("VT", 12), ("VT", 17), ("VT", 21),
                ("VT", 25), ("NM", 4), ("NM", 12), ("NM", 21), ("NM", 25), ("NM", 27),
                ("OH", 12), ("OH", 17), ("OH", 21), ("OH", 26)]
        for (iso, T) in same
            f = maximum_likelihood(SheetsGrowthModel(iso, T; data))
            @test f.estimates.r ≈ table(iso, T).r rtol = 1e-3
            @test f.estimates.K ≈ table(iso, T).K rtol = 2e-3
        end
        for (iso, T) in [("LA", 17), ("NM", 17), ("NM", 26), ("OH", 27)]   # within 2%
            f = maximum_likelihood(SheetsGrowthModel(iso, T; data))
            @test f.estimates.r ≈ table(iso, T).r rtol = 0.01
            @test f.estimates.K ≈ table(iso, T).K rtol = 0.02
        end
        # not reproduced, by us or by the authors' own code on the deposited data; the
        # values on the right are those of the R run
        R = Dict(("TN", 4) => (0.18323, 0.022500), ("OH", 4) => (0.19258, 0.064902),
                 ("VT", 4) => (0.17447, 0.54137), ("TN", 26) => (1.2372, 0.15593),
                 ("OH", 25) => (1.5215, 0.39027))
        for ((iso, T), (r, K)) in R
            f = maximum_likelihood(SheetsGrowthModel(iso, T; data))
            @test f.estimates.r ≈ r rtol = 2e-3
            @test f.estimates.K ≈ K rtol = 5e-3
            @test !isapprox(f.estimates.r, table(iso, T).r; rtol = 0.01) ||
                  !isapprox(f.estimates.K, table(iso, T).K; rtol = 0.01)
        end
        # unidentified K (exponential phase only): LA 26 °C prints K = 107300
        f = maximum_likelihood(SheetsGrowthModel("LA", 26; data))
        @test f.estimates.r ≈ table("LA", 26).r rtol = 1e-3
        @test f.estimates.K > 1e4
        # VT 26 °C: a ridge in (x0, r). The R run stops at r = 8.58 (Fig. 4 legend); our
        # optimiser goes further along the ridge to a slightly higher likelihood
        f = maximum_likelihood(SheetsGrowthModel("VT", 26; data))
        @test f.loglik >= 137.3195 - 1e-3                      # R: NLL = −137.3195
        @test f.estimates.r > 8.5
    end

    @testset "thermal responses" begin
        for iso in C.SHEETS_ISOLATES
            c = sheets_tpc(iso)
            Ts, rs = thermal_data(:sheets2021, :r; set_prefix = iso * "_")
            @test c.(Ts) ≈ rs
            @test c(19.0) ≈ (rs[3] + rs[4]) / 2
            @test c(30.0) == 0
            @test c.isolate == iso && c.source === :sheets2021
        end
        @test sheets_tpc("NM"; trait = :K)(21.0) == table("NM", 21).K
        @test sheets_tpc("OH"; trait = :zoospore_density)(4.0) == 70.3
        # "All isolates exhibited maximum viability at 21 °C" (Results); Table 2's means
        # agree except for LA, whose largest mean is at 17 °C (1.19 against 1.12)
        for iso in C.SHEETS_ISOLATES
            Ts, ms = thermal_data(:sheets2021, :MTT_mean; set_prefix = iso * "_")
            @test Ts[argmax(ms)] == (iso == "LA" ? 17.0 : 21.0)
        end
        @test_throws ArgumentError sheets_tpc("CA")
        @test_throws ArgumentError sheets_tpc("LA"; trait = :x0)
    end
    @testset "Supplementary Material: fits and profile intervals" begin
        supp(iso, T) = parameter_set(:sheets2021; set = "$(iso)_$(Int(T))C")
        # the supplement's fits are the authors' output for the deposited data; all but the
        # unidentified ones (LA 26/27, VT 26/27) are reproduced
        for iso in C.SHEETS_ISOLATES, T in C.SHEETS_TEMPERATURES
            (iso, T) in (("LA", 26.0), ("LA", 27.0), ("VT", 26.0), ("VT", 27.0)) && continue
            f = C._sheets_reference_fit(iso, T)
            @test f.estimates.r ≈ supp(iso, T).r_fit rtol = 5e-3
            @test f.estimates.K ≈ supp(iso, T).K_fit rtol = 5e-3
        end
        # profile intervals: r bounds to 1 %, K bounds to 4 % (grid interpolation in the R code)
        raw = C.TOML.parsefile(joinpath(C.COMPENDIUM_DIR, "sheets2021.toml"))["parameter"]
        interval(iso, T, sym) = only(r["uncertainty"] for r in raw
            if get(r, "set", "") == "$(iso)_$(T)C" && r["symbol"] == sym)
        for (iso, T) in (("LA", 4), ("NM", 12), ("VT", 12), ("OH", 21), ("TN", 17))
            ir = C._sheets_reference_interval(iso, T, :r)
            u = interval(iso, T, "r_fit")
            @test ir.found_lower && ir.found_upper
            @test ir.lower ≈ u["lo"] rtol = 0.01
            @test ir.upper ≈ u["hi"] rtol = 0.01
            iK = C._sheets_reference_interval(iso, T, :K)
            u = interval(iso, T, "K_fit")
            @test iK.lower ≈ u["lo"] rtol = 0.04
            @test iK.upper ≈ u["hi"] rtol = 0.04
        end
        @test !C._sheets_reference_interval("LA", 21, :r).found_upper    # supplement: FALSE
        @test !C._sheets_reference_interval("VT", 4, :K).found_upper     # supplement: FALSE
        @test_throws ArgumentError sheets_profile_interval(C._sheets_reference_fit("LA", 4), :b2)
    end
end
