module VibratoMonteCarlo
using Requires # for conditional dependencies
function __init__()
    @require DualNumbers = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74" include("deps/dual_dependencies.jl")
    @require HyperDualNumbers = "50ceba7f-c3ee-5a84-a6e8-3ad40456ec97" include("deps/hyper_dependencies.jl")
    @require ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210" include("deps/forwarddiff_dependencies.jl")
end
using FinancialMonteCarlo,Statistics;
include("utils.jl")
include("metrics/lrm.jl")
include("metrics/vibrato.jl")
include("metrics/vibrato_saltando.jl")
include("simulation/simulate_with_jump_times.jl")
include("simulation/simulate_with_last_vol.jl")
include("simulation/simulate_with_sub.jl")


export vibrato, vibrato_saltando, VibratoMonteCarloAnalytic, VibratoMonteCarloStandard

end # module

