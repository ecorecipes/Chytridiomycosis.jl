# Re-run every tutorial and report any whose rendered Markdown no longer matches the committed
# version, i.e. whose code now gives different output. CI checks that docs/src/tutorials matches
# the rendered vignettes, but nothing else re-executes them.
#
#   julia scripts/check_vignettes.jl                 # all tutorials
#   julia scripts/check_vignettes.jl 01_bsal 05_field_data
#
# Needs `quarto` and Graphviz `dot` on the PATH and the vignettes environment instantiated
# (julia --project=vignettes -e 'using Pkg; Pkg.instantiate()'). Rendering one file always
# executes its code (freeze applies only to project renders). Renders to Markdown only, in
# place, one tutorial at a time; about 15 minutes for all of them. Exits 1 if any render fails
# or any tutorial's Markdown text changed. Figures are not compared, since SVG output can differ
# byte for byte without any change in content: for tutorials whose text matched, the
# re-rendered figures and freeze caches are restored. Changed tutorials are left in the working
# tree for inspection (`git diff vignettes`); `git checkout vignettes` restores them.
# On a different platform from the one that rendered the committed files, the last printed
# digit of some floating-point results may differ: read the diff before acting on it.
const ROOT = normpath(joinpath(@__DIR__, ".."))
const VIG = joinpath(ROOT, "vignettes")

# Package precompilation messages that a first run in a fresh environment prints into a
# cell's output; they say nothing about the tutorial's results.
const NOISE = r"^\s*(Precompiling packages|\d+(\.\d+)? ms\s+✓|\d+ dependenc(y|ies) successfully precompiled)"

normalise(text) = replace(join(filter(l -> !occursin(NOISE, l), split(text, "\n")), "\n"), r"\n(\s*\n)+" => "\n\n")

# true if the rendered file differs from the committed one, ignoring precompilation messages
function differs(md)
    committed = read(Cmd(`git show HEAD:$md`; dir = ROOT), String)
    return normalise(committed) != normalise(read(joinpath(ROOT, md), String))
end

function main(names)
    all = sort(filter(d -> occursin(r"^\d\d_", d) && isfile(joinpath(VIG, d, d * ".qmd")), readdir(VIG)))
    todo = isempty(names) ? all : names
    bad = setdiff(todo, all)
    isempty(bad) || error("unknown tutorials: $(join(bad, ", "))")
    failed = String[]
    for v in todo
        t = time()
        ok = success(pipeline(Cmd(`quarto render $v/$v.qmd --to gfm`; dir = VIG);
                              stdout = devnull, stderr = devnull))
        println(rpad(v, 24), ok ? "rendered" : "RENDER FAILED", "  ", round(Int, time() - t), " s")
        ok || push!(failed, v)
    end
    mds = ["vignettes/$v/$v.md" for v in todo]
    changed = filter(md -> differs(md), mds)
    # tutorials whose text matched: restore their re-rendered figures and freeze caches, which
    # differ only byte for byte
    for v in setdiff(todo, failed)
        "vignettes/$v/$v.md" in changed && continue
        paths = filter(p -> ispath(joinpath(ROOT, p)), ["vignettes/$v/$(v)_files", "vignettes/_freeze/$v"])
        isempty(paths) || run(Cmd(`git checkout -q -- $paths`; dir = ROOT))
    end
    if isempty(changed)
        println("\nAll rendered tutorials match the committed Markdown.")
    else
        println("\nOutput changed (see `git diff vignettes`):")
        foreach(f -> println("  ", f), changed)
        run(Cmd(`git --no-pager diff --stat -- $mds`; dir = ROOT))
    end
    exit(isempty(failed) && isempty(changed) ? 0 : 1)
end

main(ARGS)
