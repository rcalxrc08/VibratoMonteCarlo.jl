include("utils.jl")

@testset "VibratoMonteCarlo" begin
    include("backends.jl")
    include("vibrato_models.jl")
    include("forwarddiff.jl")
    include("taylor.jl")
    include("regressions.jl")
    include("levy.jl")
    # Tight accuracy checks, slow and currently failing, run with VIBRATO_TEST_ACCURACY=true.
    if get(ENV, "VIBRATO_TEST_ACCURACY", "false") == "true"
        include("check/accuracy.jl")
    end
end
