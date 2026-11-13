```@raw html
<div class="fc-hero">
<div class="fc-hero-text">
<p class="fc-eyebrow"><img src="assets/logo.png" alt="" width="28" height="28"> VibratoMonteCarlo.jl</p>
```

# Smooth Monte Carlo Greeks of Any Order

```@raw html
<p class="fc-lede">Likelihood ratio, Vibrato and Vibrato Saltando estimators that are smooth in every parameter, so automatic differentiation of a Monte Carlo price gives delta, gamma and beyond, even for binary and American options.</p>
<nav class="fc-cta">
<a class="fc-btn fc-btn-primary" href="#Quick-start">Get Started <span aria-hidden="true">→</span></a>
<a class="fc-btn" href="https://github.com/rcalxrc08/VibratoMonteCarlo.jl"><svg viewBox="0 0 16 16" aria-hidden="true"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"/></svg> View on GitHub</a>
</nav>
<pre class="fc-hero-code"><code><span class="k">S0</span>    = hyper(100.0, 1.0, 1.0, 0.0)       <span class="c"># differentiate twice</span>
<span class="k">price</span> = pricer(model, zero_rate, Vibrato(mc, VibratoMonteCarloGaussLegendre()), option)
<span class="k">gamma</span> = price.epsilon12</code></pre>
</div>
<figure class="fc-hero-figure">
<img src="assets/figures/hero-surface-dark.svg" alt="Gamma of a binary call across spot and maturity, from the likelihood ratio method, with Vibrato Monte Carlo estimates as dots">
<img src="assets/figures/hero-surface-light.svg" alt="Gamma of a binary call across spot and maturity, from the likelihood ratio method, with Vibrato Monte Carlo estimates as dots">
<figcaption><span class="fc-legend"><i class="fc-dot"></i>Vibrato Monte Carlo, 100,000 paths <i class="fc-ramp"></i>gamma of a binary call, <code>PureLRM</code> with Gauss–Legendre</span>
Computed with VibratoMonteCarlo.jl when this site was built, by automatic differentiation of the prices with hyper-dual numbers. The pathwise gamma of this option is zero everywhere.</figcaption>
</figure>
</div>
```


```@raw html
<div class="fc-stats"><div class="fc-stat is-amber"><div class="fc-stat-value">2n+1<span></span></div><div class="fc-stat-label">prices give every n-th order Black–Scholes sensitivity: the Spectral Theorem</div></div><div class="fc-stat is-teal"><div class="fc-stat-value">1e-14<span></span></div><div class="fc-stat-label">accuracy of the inner Gauss–Legendre integration, split at the strikes</div></div><div class="fc-stat "><div class="fc-stat-value">6<span> models</span></div><div class="fc-stat-label">Black–Scholes, Heston, Kou, Merton, Variance Gamma, NIG</div></div><div class="fc-stat "><div class="fc-stat-value">any<span> order</span></div><div class="fc-stat-label">Greeks from HyperDualNumbers, ForwardDiff (nested), TaylorSeries, Mooncake</div></div></div>
```

## What it does

```@raw html
<p class="fc-section-lede">A plain Monte Carlo price of a binary option is a step function of the spot: its pathwise derivative is zero. <em>Conditioning on all but the last step makes the estimator smooth</em>, and automatic differentiation does the rest.</p>
<div class="fc-grid fc-features">
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><path d="M4 16 9 13l4 5 5-9 5 4 5-6"/><path d="M4 16 9 19l4-2 5 6 5-3 5 5"/><circle cx="4" cy="16" r="2"/></svg></div>
<h3>Vibrato</h3><p>Conditions on the path up to the last step and integrates the last one, for Black–Scholes, Heston and Lévy processes.</p><p class="fc-why">Smooth estimators for discontinuous payoffs, with the variance of the pathwise method.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><path d="M6 26V14M12 26V8M18 26V17M24 26V11"/><path d="M4 27h24"/></svg></div>
<h3>Vibrato Saltando <span class="fc-tag fc-tag-amber">original</span></h3><p>Conditions on the jumps, or on the subordinator, so that the inner law is exactly Gaussian, at any number of steps.</p><p class="fc-why">Exact inner integration for Kou, Merton, Variance Gamma and NIG.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><ellipse cx="16" cy="16" rx="12" ry="6"/><ellipse cx="16" cy="16" rx="7" ry="3.2"/><circle cx="16" cy="16" r="1.6" class="hl"/></svg></div>
<h3>The Spectral Theorem <span class="fc-tag fc-tag-amber">original</span></h3><p>In Black–Scholes every n-th order sensitivity is a payoff-independent combination of 2n + 1 prices; its dual links LRM and Fourier methods.</p><p class="fc-why">Ties this package to FinancialFFT.jl.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><path d="M3 26 C10 26 11 6 16 6 S22 26 29 26"/><path d="M16 6v20"/></svg></div>
<h3>Likelihood ratio</h3><p><code>PureLRM</code>: one integral against the density at maturity, from FinancialFFT.jl for Lévy models and Heston.</p><p class="fc-why">Exact up to rounding with Gauss–Legendre.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><rect x="4" y="4" width="7" height="7" rx="1.5"/><rect x="12.5" y="4" width="7" height="7" rx="1.5"/><rect x="21" y="4" width="7" height="7" rx="1.5"/><rect x="4" y="12.5" width="7" height="7" rx="1.5"/><rect x="12.5" y="12.5" width="7" height="7" rx="1.5" class="hl"/><rect x="21" y="12.5" width="7" height="7" rx="1.5"/><path d="M16 21v7"/></svg></div>
<h3>Inner integrations</h3><p>Gauss–Legendre split at the strikes, analytic, or Monte Carlo: standard, antithetic, Sobol.</p><p class="fc-why">Trade exactness for generality per method.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><circle cx="16" cy="6" r="2.5"/><circle cx="8" cy="18" r="2.5"/><circle cx="24" cy="18" r="2.5"/><path d="M16 8.5 8 15.5M16 8.5l8 7M8 20.5v6M24 20.5v6"/></svg></div>
<h3>American options</h3><p><code>Vibrato</code> over every exercise date, with Longstaff–Schwartz.</p><p class="fc-why">Greeks of early-exercise options without freezing the exercise.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><path d="M4 26 16 6l12 20z"/><path d="M10 26 16 16l6 10"/></svg></div>
<h3>Any AD backend</h3><p>DualNumbers, HyperDualNumbers, ForwardDiff (nested), ReverseDiff, TaylorSeries, Mooncake.</p><p class="fc-why">Delta, gamma and higher orders from the same price call.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><rect x="7" y="7" width="18" height="18" rx="3"/><rect x="12" y="12" width="8" height="8" rx="1"/><path d="M11 3v4M16 3v4M21 3v4M11 25v4M16 25v4M21 25v4M3 11h4M3 16h4M3 21h4M25 11h4M25 16h4M25 21h4"/></svg></div>
<h3>GPU</h3><p>The estimators run on CUDA.jl.</p><p class="fc-why">The outer simulation on the device.</p></article>
</div>
```

## How it fits together

```@raw html
<p class="fc-section-lede">The method wraps the Monte Carlo configuration of FinancialMonteCarlo.jl, so moving from plain Monte Carlo to Vibrato changes one argument.</p>
<div class="fc-arch">
<ol class="fc-flow" style="--fc-steps: 6">
<li><span class="fc-node-k">01</span><b>Model</b><span><code>BlackScholesProcess</code>, <code>KouProcess</code>, …</span></li>
<li><span class="fc-node-k">02</span><b>Outer simulation</b><span>paths up to the last step<br><code>MonteCarloConfiguration</code></span></li>
<li class="fc-core"><span class="fc-node-k">03</span><b>Conditioning</b><span><code>Vibrato</code>, <code>VibratoSaltando</code>, <code>PureLRM</code></span></li>
<li><span class="fc-node-k">04</span><b>Inner integration</b><span>Gauss–Legendre, analytic, Monte Carlo</span></li>
<li><span class="fc-node-k">05</span><b>Smooth estimator</b><span>a price that is C<sup>∞</sup> in the inputs</span></li>
<li class="fc-out"><span class="fc-node-k">06</span><b>Greeks</b><span>any order, by AD</span></li>
</ol>
</div>
```

## Results

```@raw html
<p class="fc-section-lede">Computed with the package when this site is built.</p>
<div class="fc-tabs">
<input type="radio" name="fc-results" id="fc-results-1" checked><label for="fc-results-1">Binary Greeks</label>
<input type="radio" name="fc-results" id="fc-results-2"><label for="fc-results-2">Quadrature</label>
<input type="radio" name="fc-results" id="fc-results-3"><label for="fc-results-3">Variance Gamma</label>
<div class="fc-tab-panels">
<figure class="fc-tab-panel"><img src="assets/figures/binary_greeks-dark.svg" alt="Delta and gamma of a binary call with Vibrato against Black and Scholes"><img src="assets/figures/binary_greeks-light.svg" alt="Delta and gamma of a binary call with Vibrato against Black and Scholes"><figcaption><b>Binary Greeks.</b> Delta and gamma of a binary call by Vibrato, against the closed form.</figcaption></figure>
<figure class="fc-tab-panel"><img src="assets/figures/quadrature-dark.svg" alt="Error of the inner integration"><img src="assets/figures/quadrature-light.svg" alt="Error of the inner integration"><figcaption><b>Inner integration.</b> The error of Gauss–Legendre against the other integrations.</figcaption></figure>
<figure class="fc-tab-panel"><img src="assets/figures/variance_gamma-dark.svg" alt="Variance Gamma Greeks by Vibrato Saltando"><img src="assets/figures/variance_gamma-light.svg" alt="Variance Gamma Greeks by Vibrato Saltando"><figcaption><b>Variance Gamma.</b> Vibrato Saltando on an infinite activity process.</figcaption></figure>
</div>
</div>
```

## Original research

```@raw html
<p class="fc-section-lede">Two results in this package are original work of its author, Nicola Scaramuzzino.</p>
<div class="fc-grid fc-grid-3">
<a class="fc-card fc-card-link" href="vibrato/#Vibrato-Saltando-Montecarlo"><span class="fc-card-k is-amber">Method · Master Thesis, 2017</span><h3>Vibrato Saltando</h3><p>Condition on the time of the last jump, not on a fixed step: after it the process is a Black–Scholes diffusion, and the inner expectation is known in closed form, for any number of time steps. Extended here to Variance Gamma and NIG by conditioning on the subordinator.</p></a>
<a class="fc-card fc-card-link" href="spectral/"><span class="fc-card-k is-amber">Theorem · 2024</span><h3>The Spectral Theorem</h3><p>In Black and Scholes every sensitivity of order <i>n</i> of a European option is a linear combination of the prices of 2<i>n</i> + 1 contracts, with coefficients independent of the payoff. Its dual in frequency space holds for every model with a characteristic exponent.</p></a>
</div>
```

## Quick start

```@raw html
<p class="fc-section-lede">Four steps, every output computed when this site was built.</p>
```

```@raw html
<div class="fc-step"><span class="fc-step-n">1</span><div><h3>Install the package</h3><p>Some dependencies are in a registry of their own.</p></div></div>
```

```julia
import Pkg
Pkg.Registry.add(Pkg.RegistrySpec(url = "https://github.com/rcalxrc08/LocalRegistry.git"))
Pkg.Registry.add("General")
Pkg.add("VibratoMonteCarlo")
```

```@raw html
<div class="fc-step"><span class="fc-step-n">2</span><div><h3>Mark the input to differentiate</h3><p>A hyper-dual spot carries delta and gamma through the price.</p></div></div>
```

```@example quickstart
using VibratoMonteCarlo, FinancialMonteCarlo, FinancialToolbox, HyperDualNumbers

S0 = hyper(100.0, 1.0, 1.0, 0.0)
K, r, d, T, σ = 100.0, 0.02, 0.01, 1.0, 0.2
model = BlackScholesProcess(σ, Underlying(S0, d))
zero_rate = ZeroRate(r)
mc = MonteCarloConfiguration(100_000, 10; sampling = FinancialMonteCarlo.AntitheticMC())
nothing # hide
```

```@raw html
<div class="fc-step"><span class="fc-step-n">3</span><div><h3>Price with Vibrato</h3><p>Value, delta and gamma from one evaluation, against the formula of Black and Scholes.</p></div></div>
```

```@example quickstart
price = pricer(model, zero_rate, Vibrato(mc, VibratoMonteCarloGaussLegendre()), EuropeanOption(T, K))
exact = blsprice(S0, K, r, T, σ, d)
(vibrato = (price.value, price.epsilon1, price.epsilon12), exact = (exact.value, exact.epsilon1, exact.epsilon12))
```

```@raw html
<div class="fc-step"><span class="fc-step-n">4</span><div><h3>A binary option on a jump model</h3><p>Where the pathwise gamma is zero, with Vibrato Saltando.</p></div></div>
```

```@example quickstart
kou = KouProcess(0.2, 1.0, 0.3, 20.0, 10.0, Underlying(S0, d))
price = pricer(kou, zero_rate, VibratoSaltando(mc, VibratoMonteCarloGaussLegendre()), BinaryEuropeanOption(T, K))
(value = price.value, delta = price.epsilon1, gamma = price.epsilon12)
```

```@raw html
<div class="fc-next">
<a href="methods/"><b>Methods</b><span>The estimators, their inner integrations, how to choose.</span></a>
<a href="vibrato/"><b>Vibrato and Vibrato Saltando</b><span>Conditioning, with the original Saltando method.</span></a>
<a href="spectral/"><b>The Spectral Theorem</b><span>2n + 1 prices for every n-th order Greek.</span></a>
<a href="api/"><b>API reference</b><span>Every exported name.</span></a>
</div>
```

## The ecosystem

The package adds methods to the `pricer` of [FinancialMonteCarlo.jl](https://github.com/rcalxrc08/FinancialMonteCarlo.jl),
as [FinancialFFT.jl](https://github.com/rcalxrc08/FinancialFFT.jl) does, and takes the densities of the Lévy processes
from FinancialFFT.jl. The implementation is based on the Master Thesis of Nicola Scaramuzzino [scaramuzzino2017](@cite)
and its subsequent development. Tests: `] test VibratoMonteCarlo`; the GPU tests run with
`julia --project=test/cuda -e 'import Pkg; Pkg.instantiate(); include("test/cuda/runtests.jl")'`. More projects at
[rcalxrc08.github.io](https://rcalxrc08.github.io/).
