module VibratoMonteCarlo
using Requires # for conditional dependencies
function __init__()
    @require DualNumbers = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74" include("deps/dual_dependencies.jl")
    @require HyperDualNumbers = "50ceba7f-c3ee-5a84-a6e8-3ad40456ec97" include("deps/hyper_dependencies.jl")
    @require ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210" include("deps/forwarddiff_dependencies.jl")
    @require ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267" include("deps/reversediff_dependencies.jl")
    @require TaylorSeries = "6aa5eb33-94cf-58f4-a9d0-e4b2c4fc25ea" include("deps/taylorseries_dependencies.jl")
    @require FastGaussQuadrature = "442a2c76-b920-505d-bb47-c5924d526838" include("deps/fastgaussquadrature.jl")
end
using FinancialMonteCarlo, Statistics, MuladdMacro, IrrationalConstants
include("utils.jl")
include("metrics/lrm.jl")
include("metrics/lrm_aug.jl")
include("metrics/vibrato.jl")
include("metrics/pure_lrm.jl")
include("metrics/vibrato_saltando.jl")
include("simulation/simulate_with_jump_times.jl")
include("simulation/simulate_with_last_vol.jl")
include("simulation/simulate_with_sub.jl")

export vibrato, vibrato_saltando, pure_lrm, VibratoMonteCarloAnalytic, VibratoMonteCarloStandard

end # module
