# Pricing methods for FinancialMonteCarlo's interface: pricer(mcProcess, rfCurve, method, payoff), as FinancialFFT's LewisMethod.
import FinancialMonteCarlo: pricer
abstract type AbstractVibratoMethod <: FinancialMonteCarlo.AbstractMethod end
# A method is a scalar under broadcasting, like AbstractVibrato and FinancialFFT's methods: pricer.(models, rfCurve, method, payoffs)
# prices every model and payoff with the same method instead of trying to iterate over it.
Base.broadcastable(x::AbstractVibratoMethod) = Ref(x)

# Payoffs that depend on S_T only: the last step (or the whole maturity) is integrated against its density.
const EuropeanLike = Union{FinancialMonteCarlo.EuropeanPayoff, FinancialMonteCarlo.AbstractSmile}

"""
    Vibrato(mc, vb_mc)

Vibrato Monte Carlo: `mc` (a `MonteCarloConfiguration`) simulates the path up to the last time step, and `vb_mc`
(`VibratoMonteCarloGaussLegendre`, `VibratoMonteCarloAnalytic`, `VibratoMonteCarloStandard`) integrates the payoff
against the density of the last step with the likelihood ratio method.

`pricer(mcProcess, rfCurve, Vibrato(mc, vb_mc), payoff)`:
- `BlackScholesProcess` and `HestonProcess`: Gaussian last step (for Heston with the simulated variance).
- Lévy processes: FFT density of the last step (requires `VibratoMonteCarloAnalytic`).
- an `AmericanPayoff` (`AmericanOption`, `BinaryAmericanOption`): exercise at the `Nstep` grid dates, policy from
  Longstaff-Schwartz. On `BlackScholesProcess` vibrato over every exercise date; on Lévy processes and
  `HestonProcess` the exercise decisions are frozen when differentiating, so second order sensitivities are biased.
"""
struct Vibrato{M <: FinancialMonteCarlo.AbstractMonteCarloConfiguration, V <: AbstractVibrato} <: AbstractVibratoMethod
    mc::M
    vb_mc::V
end

"""
    VibratoSaltando(mc, vb_mc)

Vibrato conditioned on the jumps: on the last jump time for finite activity processes (`KouProcess`, `MertonProcess`),
on the subordinator `G_T` for `VarianceGammaProcess` and `NormalInverseGaussianProcess`. The law of `log S_T` given that
is Gaussian, and `vb_mc` integrates the payoff against it. Exact for any `Nstep`. On VG and NIG it also differentiates over
`κ` and `T`; on Kou and Merton it does not support derivatives over `λ` or `T` (it warns).
"""
struct VibratoSaltando{M <: FinancialMonteCarlo.AbstractMonteCarloConfiguration, V <: AbstractVibrato} <: AbstractVibratoMethod
    mc::M
    vb_mc::V
end

"""
    PureLRM(mc, vb_mc)

Likelihood ratio method on the density at maturity, with no path simulation: Gaussian for Black-Scholes, from the FFT
for the other processes (requires `VibratoMonteCarloAnalytic`). `mc` only provides the seed and, for
`VibratoMonteCarloStandard`, the random numbers.
"""
struct PureLRM{M <: FinancialMonteCarlo.AbstractMonteCarloConfiguration, V <: AbstractVibrato} <: AbstractVibratoMethod
    mc::M
    vb_mc::V
end

# Not exported, mainly test oracles: the same paths as Vibrato (Black-Scholes) and VibratoSaltando, with the conditional
# expectation in closed form (European and binary European options only) and pathwise derivatives.
"""
    VibratoMonteCarlo.Conditional(mc)

Same paths as `Vibrato` on `BlackScholesProcess`, with the last step priced in closed form (`blsprice`, `blsbin`).
"""
struct Conditional{M <: FinancialMonteCarlo.AbstractMonteCarloConfiguration} <: AbstractVibratoMethod
    mc::M
end

"""
    VibratoMonteCarlo.ConditionalSaltando(mc)

Same paths as `VibratoSaltando`, with what the conditioning leaves priced in closed form: the exact inner integration limit
of `VibratoSaltando`.
"""
struct ConditionalSaltando{M <: FinancialMonteCarlo.AbstractMonteCarloConfiguration} <: AbstractVibratoMethod
    mc::M
end
