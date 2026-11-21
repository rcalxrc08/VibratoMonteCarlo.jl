module VibratoMonteCarlo
# AD backends are supported through package extensions, see ext/. The stages, in order (contracts in src/0N_*/CONTEXT.md):
# methods and inner integrations, outer simulation, inner integration of the last step, European pricers, American pricers.
using FinancialMonteCarlo, Statistics, MuladdMacro, IrrationalConstants, PatternVectors
using FastGaussQuadrature: gausshermite, gausslegendre
include("01_methods/utils.jl")
include("01_methods/methods.jl")
include("02_simulation/simulate_with_jump_times.jl")
include("02_simulation/simulate_with_last_vol.jl")
include("02_simulation/subordinator.jl")
include("03_inner_integration/lrm.jl")
include("03_inner_integration/lrm_aug.jl")
include("04_european/vibrato.jl")
include("04_european/pure_lrm.jl")
include("04_european/vibrato_saltando.jl")
include("05_american/vibrato_am.jl")

export pricer, Vibrato, VibratoSaltando, PureLRM, VibratoMonteCarloAnalytic, VibratoMonteCarloStandard, VibratoMonteCarloGaussLegendre

end # module
