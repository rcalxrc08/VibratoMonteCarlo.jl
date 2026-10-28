module VibratoMonteCarlo
# AD backends and FastGaussQuadrature are supported through package extensions, see ext/.
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
