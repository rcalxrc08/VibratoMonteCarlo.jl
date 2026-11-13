```@raw html
<div class="fc-page-head">
<p class="fc-kicker"><span>Theory 04</span> Variance and accuracy</p>
```

# Numerical refinements

The estimators of [Conditioning: Vibrato and Vibrato Saltando](@ref) are exact in theory; the following refinements make the inner expectation accurate in practice.
In all of them the inner law is Gaussian, ``X_T | X_{\tau} \sim \mathcal{N}(\mu, s^2)``, and the inner expectation is computed on the standardized variable ``Z = (X_T - \mu)/s``.

```@raw html
<div class="fc-meta"><span>centering</span><span>inner samples</span><span>quadrature</span></div>
<nav class="fc-toc">
<a href="#Payoff-Centering"><b>01</b>Payoff Centering</a>
<a href="#Independent-Inner-Samples"><b>02</b>Independent Inner Samples</a>
<a href="#Quadrature"><b>03</b>Quadrature</a>
</nav>
</div>
```

## Payoff Centering
For every ``\theta`` the density integrates to one:

```math
\int_{\mathbb{R}} f(x,\theta) dx = 1
```

hence all of its derivatives vanish and, for ``n \geq 1``:

```math
\mathbb{E}(B_n(\{\partial^{i}_{\theta} \log{f(X_T,\theta)}\}_{i = 1 : n})) = 0
```

Therefore any constant ``c`` not depending on ``\theta`` can be subtracted from the payoff without changing the sensitivities:

```math
\partial^{n}_{\theta} V_0 = \mathbb{E}((\Phi(X_T) - c) B_n(\{\partial^{i}_{\theta} \log{f(X_T,\theta)}\}_{i = 1 : n}))
```

We take ``c = \Phi(\mu)``, the payoff at the centre of the inner density, and add it back to the price.
The price is unchanged, while the score terms are now multiplied by ``\Phi(\mu + s Z) - \Phi(\mu) \approx \Phi'(\mu) s Z`` instead of ``\Phi(\mu + s Z)``.
In vibrato ``s = \sigma \sqrt{T - \tau}`` is small and the score terms grow like ``s^{-n}``, so without centering the variance of the delta estimator grows like ``s^{-2}``; with centering it stays bounded.
Together with antithetic pairs ``\pm Z``, the second order term becomes ``\frac{1}{2}(\Phi(\mu + s Z) + \Phi(\mu - s Z)) - \Phi(\mu) \approx \frac{1}{2}\Phi''(\mu) s^2 Z^2``, which compensates the ``s^{-2}`` of the gamma score.
In practice this reduces the noise of the gamma estimated with inner Monte Carlo by a factor 10 to 80.

## Independent Inner Samples
The inner expectation is estimated with ``M`` samples ``Z_1, \dots, Z_M`` for each of the ``N`` outer paths.
If the same samples were used for every outer path, the inner error would be the same for all paths and would not average out: the estimator would converge, as ``N \to \infty``, to a biased value with error of order ``M^{-1/2}`` (deltas 8-15% off with ``M = 200``).
Drawing new inner samples for every outer path makes the inner errors independent, and the total error decreases with both ``N`` and ``M``.
With quasi Monte Carlo inner samples, the same Sobol points are randomly shifted modulo 1 for every outer path (Cranley-Patterson rotation [cranley1976](@cite)), which keeps each inner estimate unbiased.

## Quadrature
Alternatively, the inner expectation can be computed as a deterministic integral on a uniform grid in ``Z``.
The grid must reach far enough in the tails: for a call the payoff grows exponentially in the upper tail, and truncating at ``\pm 5`` standard deviations loses enough of the second order score terms to bias the volga by 3% to 45% on Variance Gamma. The default grid spans ``\pm 7`` standard deviations.
