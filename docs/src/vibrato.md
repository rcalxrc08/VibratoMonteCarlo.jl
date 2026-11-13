```@raw html
<div class="fc-page-head">
<p class="fc-kicker"><span>Theory 03</span> Conditioning, and the original Vibrato Saltando</p>
```

# Conditioning: Vibrato and Vibrato Saltando


```@raw html
<div class="fc-meta"><span>Vibrato</span><span>Vibrato Saltando</span><span>Lévy</span></div>
<nav class="fc-toc">
<a href="#Conditional-Expectation"><b>01</b>Conditional Expectation</a>
<a href="#Conditional-Saltando-Expectation"><b>02</b>Conditional Saltando Expectation</a>
<a href="#Vibrato-Montecarlo-Method"><b>03</b>Vibrato Montecarlo Method</a>
<a href="#Vibrato-Saltando-Montecarlo"><b>04</b>Vibrato Saltando Montecarlo</a>
</nav>
</div>
```
## Conditional Expectation
The payoff is not smooth, but its conditional expectation is. Let's fix ``\tau < T``, let's call ``h = T - \tau`` and let's assume that the underlying process is a Levy process. Since:

```math
V_0=\mathbb{E}(\Phi(X_T))=\mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))
```

```math
\partial^{n}_{\theta} V_0=\partial^{n}_{\theta} \mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))=\mathbb{E}(\partial^{n}_{\theta} \mathbb{E}(\Phi(X_{T})|X_{\tau}))
```

The inner expectation is a smooth function of ``X_{\tau}`` and of the parameters, even if ``\Phi`` is a kink or a jump. Hence, when it is known in closed form, automatic differentiation can go through it path by path, at any order.

### Black and Scholes
The inner expectation is the price of the same option over the remaining time ``h``, for a call:

```math
C(S,h) = S e^{-d h} N(d_1) - K e^{-r h} N(d_2)
```

```math
d_1 = \frac{\log(\frac{S}{K}) + (r-d+\frac{\sigma^2}{2}) h}{\sigma \sqrt{h}}, \qquad d_2 = d_1 - \sigma \sqrt{h}
```

where ``N`` is the standard normal distribution function and ``\varphi`` its density. Discounting up to ``\tau``, the method `Conditional` computes:

```math
V_0 = e^{-r \tau}\, \mathbb{E}(C(S_{\tau},h))
```

The derivatives of the inner price with respect to ``S``, generated with Symbolics, for a call and for a binary call:
```@eval
Main.DocumentationFormulas.conditional_derivatives()
```
Since ``\partial_{S_0} S_{\tau} = \frac{S_{\tau}}{S_0}``, the delta is then:

```math
\Delta = e^{-r \tau}\, \mathbb{E}(\frac{S_{\tau}}{S_0} \partial_{S} C(S_{\tau},h))
```

This is the best one can do when it applies, but it needs a closed form for the inner expectation, so it is limited to European and binary options on Gaussian steps. Vibrato removes this limit by computing the same inner expectation by integration, with no closed form.

## Conditional Saltando Expectation
Let's assume now that ``X_t`` is a finite activity Levy process, like Kou [kou2002](@cite) or Merton [merton1976](@cite), with jump intensity ``\lambda`` and ``\kappa = \mathbb{E}(e^{Y_1} - 1)``, where ``Y_i`` are the jump sizes. Instead of a fixed ``\tau``, let's take the time of the last jump before ``T``, with ``\tau = 0`` if there are no jumps.

This ``\tau`` is not a stopping time, since knowing it needs a look at the future of the path, so we condition on all the jumps and on the path up to ``\tau``. After ``\tau`` there are no more jumps and the process is a Brownian motion with a modified drift, which means a Black and Scholes model with dividend yield ``d + \lambda \kappa``. Hence the inner expectation is known in closed form, even if the process has jumps:

```math
V_0 = \mathbb{E}(e^{-r \tau} C_{\lambda \kappa}(S_{\tau},T - \tau))
```

where ``C_{\lambda \kappa}`` is the price of the previous section with ``d`` replaced by ``d + \lambda \kappa``. The discount is inside the expectation since ``\tau`` is different for every path. This is what `ConditionalSaltando` computes.
With ``h = T - \tau``, the same derivatives become:
```@eval
Main.DocumentationFormulas.conditional_saltando_derivatives()
```
On the same paths, this is the exact limit of Vibrato Saltando for a perfect inner integration. It has the same limits of the fixed ``\tau`` case, it works only for European and binary options, and the same ones of Vibrato Saltando: ``\lambda`` and ``T`` can't be differentiated, since the last jump time depends on them.

For Variance Gamma and Normal Inverse Gaussian, `ConditionalSaltando` conditions on the subordinator ``G_T`` instead (see Vibrato Saltando below): the inner price is Black and Scholes with volatility ``\sigma \sqrt{G_T / T}`` and the rest of the drift as a dividend yield. There the price is a one dimensional integral over ``G_T``, the formula of [madan1998](@cite) for Variance Gamma, so the estimator is mostly a reference for `VibratoSaltando`.

## Vibrato Montecarlo Method
The pathwise method needs a smooth payoff, while the LRM has weights whose variance grows with the order of the sensitivity. Vibrato uses the best of each: it differentiates the path pathwise up to a time ``\tau < T``, and applies the LRM only to what is left. The method was introduced in [giles2009](@cite), and its combination with automatic differentiation for high order sensitivities is studied in [pages2018](@cite).

Let's fix ``\tau < T``, let's call ``h = T - \tau`` and let's assume that the underlying process is a Levy process. Since:

```math
V_0=\mathbb{E}(\Phi(X_T))=\mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))
```

```math
\partial^{n}_{\theta} V_0=\mathbb{E}(\partial^{n}_{\theta} \mathbb{E}(\Phi(X_T)|X_{\tau}))
```

The inner expectation is smooth in ``\theta`` even if ``\Phi`` is not, since it is an integral of ``\Phi`` against a smooth density. This is why the outer expectation can be differentiated. Applying the LRM to the inner expectation:

```math
\partial^{n}_{\theta} V_0=\mathbb{E}(\mathbb{E}(\Phi(X_{T}) B_n(\{\partial^{i}_{\theta} \log{f_{X_T|X_{\tau}}(X_T,\theta)}\}_{i = 1 : n})|X_{\tau}))
```

As the traditional LRM, this needs the transition density. It is known for Ito processes, where the last step is Gaussian. For more general Levy processes it has to be inverted numerically from the characteristic function.

### Gaussian last step
If the last step is Gaussian, ``X_T|X_{\tau} \sim \mathcal{N}(\mu, s^2)``, and ``\mu`` and ``s`` carry all the dependence on the parameters, with:

```math
Z = \frac{X_T - \mu}{s}
```

then the scores with respect to ``\mu`` and ``s``, generated from `log_density`, are:
```@eval
Main.DocumentationFormulas.gaussian_scores()
```
and the score with respect to a generic parameter follows from the chain rule:

```math
\partial_{\theta} \log{f} = \partial_{\mu} \log{f}\, \partial_{\theta} \mu + \partial_{s} \log{f}\, \partial_{\theta} s
```

In case of the Black and Scholes model:

```math
\mu = X_{\tau} + (r-d-\frac{\sigma^2}{2}) h, \qquad s = \sigma \sqrt{h}
```

Note that ``\mu`` depends on the parameters also through the path: this is the pathwise part of the method.

### Why automatic differentiation works
The weights are never written by hand in the package. Let's call ``\mathrm{v}(a)`` the value of a number ``a`` without its derivatives. Then the sample points carry no derivative:

```math
x(z) = \mathrm{v}(\mu) + \mathrm{v}(s)\, z
```

so the payoff is never differentiated, and all of the parameters enter through the density, in the factor:

```math
m(g) = e^{g - \mathrm{v}(g)}
```

which has value one and, by Faa di Bruno's formula, derivatives ``\partial^{n}_{\theta} m(g) = B_n(\{\partial^{i}_{\theta} g\}_{i = 1 : n})``. Hence evaluating at ``g = \log{f(x(z);\mu,s)}`` gives the weights of every order out of any automatic differentiation number. The estimator on ``N`` paths is, with ``\mu_j`` and ``s_j`` the values of the path ``j`` and ``c_j = \Phi(\mathrm{v}(\mu_j))`` the payoff at the centre of the density (see [Payoff Centering](@ref)):

```math
\hat{V}_0 = \frac{e^{-r T}}{N} \sum_{j=1}^{N} \int_{\mathbb{R}} \left((\Phi(x_j(z)) - c_j)\, m(\log{f(x_j(z);\mu_j,s_j)}) + c_j\right) \varphi(z) dz
```

where ``\varphi`` is the standard normal density and ``x_j(z) = \mathrm{v}(\mu_j) + \mathrm{v}(s_j) z``. The value is the plain Monte Carlo price, and its derivatives are the vibrato estimators.

### Weights for Black and Scholes
Writing the path as ``X_{\tau} = \log(S_0) + (r-d-\frac{\sigma^2}{2}) \tau + \sigma W_{\tau}`` and ``\delta = \sqrt{h}``, we generate the weights with [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), starting from the same `log_density` used by the package and evaluating them at ``X_T = \mu + \sigma \delta Z``:
```@eval
Main.DocumentationFormulas.vibrato_weights()
```
Each sensitivity is then, with ``w`` the corresponding weight:

```math
\partial^{n}_{\theta} V_0 = e^{-r T}\, \mathbb{E}(\mathbb{E}((\Phi(\mu + \sigma \delta Z) - \Phi(\mu))\, w|X_{\tau}))
```

For ``\tau = 0`` these are the LRM weights of [Likelihood Ratio Method for Black and Scholes](@ref). The smaller ``\delta`` is, the larger the weights are, and only the centering of the payoff keeps their variance under control. In ``w_{\mathcal{V}}`` the term with ``W_{\tau}`` comes from the path, while ``Z^2 - 1`` comes from the last step alone.

## Vibrato Saltando Montecarlo

!!! tip "Original contribution"
    Vibrato Saltando is an original method of the author of this package, Nicola Scaramuzzino, introduced in the
    Master Thesis [scaramuzzino2017](@cite) and extended here to infinite activity processes.
Vibrato needs the law of the last step, and for a Levy process that is only available numerically. For Variance Gamma it can even be singular on a short step.
Vibrato Saltando [scaramuzzino2017](@cite) chooses what to condition on so that the inner law is exactly Gaussian, for any number of time steps.

The method was introduced in the Master Thesis [scaramuzzino2017](@cite), Section 7.6, for finite activity Levy processes. "Vibrato" comes from the method of Giles [giles2009](@cite); "Saltando" (Italian for "jumping") refers both to the jumps (*salti*) of the processes it was made for and to a bowing technique of string players: to change note across the instrument, instead of keeping the bow on the strings they lift it and let it bounce on the next string.

### Finite activity processes
Let's assume that ``X_t`` is a finite activity Levy process, like Kou or Merton:

```math
X_t = \log(S_0) + (r-d-\frac{\sigma^2}{2}-\lambda \kappa) t + \sigma W_t + \sum_{i=1}^{N_t} Y_i, \qquad \kappa = \mathbb{E}(e^{Y_1} - 1)
```

and let's call ``\tau`` the time of the last jump before ``T``, with ``\tau = 0`` if there are no jumps.
This ``\tau`` is random and it is not a stopping time, since knowing it needs a look at the future of the path. So we condition on all the jumps and on ``W_{\tau}``: the jumps are independent of ``W``, and after ``\tau`` there is nothing left but the diffusion. Calling ``J = \sum_{i=1}^{N_T} Y_i``:

```math
X_T|\tau, J, W_{\tau} \sim \mathcal{N}(\mu, s^2)
```

```math
\mu = X_{\tau} + (r-d-\frac{\sigma^2}{2}-\lambda \kappa)(T - \tau), \qquad s = \sigma \sqrt{T - \tau}
```

The method is the one of the previous section, with a step that changes from path to path: the jumps go in the pathwise part, and the LRM is applied to the Gaussian diffusion after the last jump. Hence the payoff is integrated exactly, and ``Z`` has no discretization error. On a path with no jumps ``\tau = 0``, and Vibrato Saltando is the plain LRM over the whole maturity [scaramuzzino2017](@cite).
Generating the weights as before, with ``\delta = \sqrt{T - \tau}``:
```@eval
Main.DocumentationFormulas.saltando_weights()
```
They are the same weights of Black and Scholes, since neither ``J`` nor ``\lambda \kappa`` depend on ``S_0`` or ``\sigma``. What changes is ``\delta``, which is now different for every path.

Jumps are applied at the next point of the time grid, so the simulation returns the grid point before the step containing the last jump. This keeps the last step exact and at least one time step long: otherwise ``T - \tau`` could be zero, with infinite weights.

The last jump time depends on ``\lambda`` and the length of the last step on ``T``, while the weights above consider them fixed: differentiating with respect to ``\lambda`` or ``T`` is not supported.

### Infinite activity processes
Infinite activity processes like Variance Gamma [madan1998](@cite) and Normal Inverse Gaussian [barndorff1997](@cite) have no last jump, but they can be written as a Brownian motion with drift evaluated at a random time ``G_t`` (the subordinator):

```math
X_t = \log(S_0) + (r - d - \psi) t + \beta G_t + \sigma W_{G_t}
```

where ``\beta`` and ``\sigma`` are model parameters (called `θ` and `σ` in FinancialMonteCarlo), and ``\psi`` is the drift correction such that ``\mathbb{E}(S_T) = S_0 e^{(r-d)T}``.
Conditioning on the subordinator, the inner law is exactly Gaussian:

```math
X_T|G_T \sim \mathcal{N}(\mu, s^2)
```

```math
\mu = \log(S_0) + (r - d - \psi) T + \beta G_T, \qquad s = \sigma \sqrt{G_T}
```

Hence only ``G_T`` has to be simulated (Gamma for Variance Gamma, Inverse Gaussian for Normal Inverse Gaussian), and the LRM applies to the Gaussian inner expectation with no time discretization.
Since ``G_T`` is one dimensional, it is drawn by inversion with antithetic and Sobol sampling: the quantiles of its law at ``U`` and ``1 - U``, or at the points of a one dimensional Sobol sequence.
The delta and gamma weights are once more the same, now with ``\delta = \sqrt{G_T}``:
```@eval
Main.DocumentationFormulas.subordinated_weights()
```
For the vega, ``\psi`` depends on ``\sigma`` too: there the automatic differentiation of the factor ``m`` takes care of it.

This is the extension to more general models that the thesis left as future work [scaramuzzino2017](@cite). Unlike ``\lambda`` for finite activity processes, the parameters of the law of ``G_T`` (its variance rate ``\kappa``, and ``T``) can be differentiated: ``G_T`` is drawn from their values, and they enter through its density, in the same factor ``m`` applied to the outer expectation:

```math
V_0 = \mathbb{E}(\mathbb{E}(e^{-rT} \Phi(X_T)|G_T)\, m(\log{f_{G_T}(G_T;\kappa,T)}))
```

The draws of a Gamma or Inverse Gaussian sampler are not smooth in ``\kappa``, while the density is. The noise of these derivatives is reduced by subtracting the sample mean of the inner expectations, as in [Payoff Centering](@ref).

### Summary
All of the methods condition on something, and integrate the last part against its law. They differ in what they condition on, which gives the law of what is left (``\mu`` is its mean), and in how the inner expectation is computed:

| Process | Method | Conditioning on | Law of ``X_T`` | Inner expectation |
|:--------|:---------|:----------------|:------------------------------------|:------------------|
| Black and Scholes | `Conditional` | path up to ``T - h`` | ``\mathcal{N}(\mu, \sigma^2 h)`` | closed form |
| Black and Scholes | `Vibrato` | path up to ``T - h`` | ``\mathcal{N}(\mu, \sigma^2 h)`` | integration |
| Heston [heston1993](@cite) | `Vibrato` | path and variance ``v`` up to ``T - h`` | ``\mathcal{N}(\mu, v_{T-h} h)`` | integration |
| Levy | `Vibrato` | path up to ``T - h`` | FFT density | integration |
| Kou, Merton | `ConditionalSaltando` | jumps and path up to the last jump ``\tau`` | ``\mathcal{N}(\mu, \sigma^2 (T - \tau))`` | closed form |
| Kou, Merton | `VibratoSaltando` | jumps and path up to the last jump ``\tau`` | ``\mathcal{N}(\mu, \sigma^2 (T - \tau))`` | integration |
| Variance Gamma, Normal Inverse Gaussian | `ConditionalSaltando` | subordinator ``G_T`` | ``\mathcal{N}(\mu, \sigma^2 G_T)`` | closed form |
| Variance Gamma, Normal Inverse Gaussian | `VibratoSaltando` | subordinator ``G_T`` | ``\mathcal{N}(\mu, \sigma^2 G_T)`` | integration |
