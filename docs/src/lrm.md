# Pathwise and likelihood ratio methods

## The basic settings
Let's assume we have an underlying stock price varying as a stochastic process called ``S_t`` and
let's assume without loss of generality that ``S_t>0`` a.e..\\
In that case we can define as ``X_t`` the following:

```math
X_t=\log(S_t)
```

Let's assume we want to price and compute the sensitivities of a european option with maturity ``T`` and payoff ``\phi``.

We can express the price of such option as:

```math
V_0=\mathbb{E}(\phi(S_T))
```

or, in function of ``X_T``:

```math
V_0=\mathbb{E}(\Phi(X_T))
```

Where:

```math
\Phi(x)=\phi(e^x)
```

## Pathwise Method
The pathwise method and the likelihood ratio method are the two classic direct estimators of sensitivities [broadie1996](@cite), see also [glasserman2004](@cite).

Since:

```math
V_0=\mathbb{E}(\Phi(X_T))
```

```math
\partial_\theta V_0=\partial_\theta \mathbb{E}(\Phi(X_T)) = \mathbb{E}(\dfrac{\partial \Phi(X_T)}{\partial X_T}\dfrac{\partial X_T}{\partial \theta})
```

Where the following is called tangent process:

```math
Y_t=\dfrac{\partial X_t}{\partial \theta}
```

In case of a European Option with maturity ``T`` and strike price ``K``:

```math
\Phi(x) = max(e^x-K,0)
```

and:

```math
\dfrac{\partial \Phi(x)}{\partial x} = e^x \mathbb{1}_{ e^{x} > K } (x)
```

Hence:

```math
\partial_\theta V_0=\mathbb{E}(\dfrac{\partial \Phi(X_T)}{\partial X_T}\dfrac{\partial X_T}{\partial \theta})
```

Unfortunately this method does not provide usefull results for binary options or for n-th order sensitivities, indeed in case of binary options we have:

```math
\Phi(x)=\mathbb{1}(x)
```

Hence:

```math
\partial_x \Phi(x)=0
```

## Likelihood Ratio Method
Starting from the price of the basic settings, if we develop the integral:

```math
V_0=\mathbb{E}(\Phi(X_T))=\int_{\mathbb{R}}  \Phi(x) f(x) dx
```

where the following is the density of the log yield:

```math
f(x)=f_{X_T}(x)
```

If we differentiate for a generic parameter ``\theta``:

```math
\partial_\theta V_0=\partial_\theta \int_{\mathbb{R}}  \Phi(x) f(x,\theta) dx = \int_{\mathbb{R}}  \Phi(x) \partial_\theta f(x,\theta) dx
```

or:

```math
\partial_\theta V_0 = \int_{\mathbb{R}}  \Phi(x) \partial_\theta(\log(f(x,\theta))) f(x,\theta) dx
```

or:

```math
\partial_\theta V_0 = \mathbb{E}(\Phi(X_T) \partial_\theta \log(f(X_T,\theta)))
```

In case of n-th order sensitivities:

We observe that for a generic positive differentiable function f we can have:

```math
g(x,\theta) =  \log{f(x,\theta)}
```

or:

```math
f(x,\theta) =  e^{g(x,\theta)}
```

Then, according to [Faa di Bruno's formula](https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula):

```math
\partial^{n}_{\theta} e^{g(x,\theta)} = e^{g(x,\theta)} B_n(\{\partial^{i}_{\theta} g(x,\theta)\}_{i = 1 : n})
```

Where ``B_n`` is the complete exponential [Bell polynomial](https://en.wikipedia.org/wiki/Bell_polynomials) of n-th order.
### Bell polynomials up to order 6
```@eval
Main.DocumentationFormulas.bell_polynomials()
```

Hence:

```math
\partial^{n}_{\theta} f(x,\theta) = f(x,\theta) B_n(\{\partial^{i}_{\theta} \log{f(x,\theta)}\}_{i = 1 : n})
```

By applying this to the original problem:

```math
\partial^{n}_{\theta} V_0=\partial^{n}_{\theta} \int_{\mathbb{R}}  \Phi(x) f(x,\theta) dx = \int_{\mathbb{R}}  \Phi(x) \partial^{n}_{\theta} f(x,\theta) dx
```

By applying the previous result:

```math
\partial^{n}_{\theta} V_0= \int_{\mathbb{R}}  \Phi(x) f(x,\theta) B_n(\{\partial^{i}_{\theta} \log{f(x,\theta)}\}_{i = 1 : n}) dx
```

or, in terms of expectation:

```math
\partial^{n}_{\theta} V_0 = \mathbb{E}(\Phi(X_T) B_n(\{\partial^{i}_{\theta} \log{f(X_T,\theta)}\}_{i = 1 : n}))
```

### Likelihood Ratio Method for Black and Scholes
In case of trivial models the density is known in analytic form, otherwise one need to compute it numerically from the characteristic function.
In case of the Black and Scholes Model [blackscholes1973](@cite), and a European Option we have the following:

```math
X_t=\log(S_0)+(r-d-\dfrac{\sigma^2}{2}) t + \sigma W_t
```

```math
\log(f(x,S_0,r,d,\sigma,T)) = -\dfrac{1}{2}((\dfrac{x-(\log(S_0)+(r-d-\dfrac{\sigma^2}{2}) T)}{\sigma \sqrt{T}})^2 + \log(2 \pi)) - \log(\sigma\sqrt{T})
```

From this we can compute the various partial derivatives. Calling ``\delta = \sqrt{T}`` and:

```math
Z = \frac{x - \log(S_0) - (r-d-\frac{\sigma^2}{2}) T}{\sigma \delta}
```

which is a standard normal when ``x = X_T``, [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) generates from the `log_density` of the package:
```@eval
Main.DocumentationFormulas.lrm_scores()
```

Hence, a delta for a generic european option can be computed as the following expectation:

```@eval
Main.DocumentationFormulas.lrm_delta()
```

or, in terms of ``W_T``:

```@eval
Main.DocumentationFormulas.lrm_delta_wt()
```

Let's assume now that ``\Phi`` is a forward i.e.:

```math
\Phi(x)=e^x
```

Since ``\Delta = e^{(r-d)T}``:

```math
\mathbb{E}(e^{-\frac{\sigma^2}{2} T + \sigma W_T}\,W_T) = \sigma T
```

### Spectral Theorem for Black and Scholes
The scores above are polynomials of degree at most 2 in ``x``, so the weight of order ``n`` is a polynomial of degree
``2n`` in ``X_T``, and every sensitivity of order ``n`` is a combination, with coefficients that do not depend on the
payoff, of the prices of the ``2n + 1`` contracts ``\mathbb{E}(\Phi(X_T) X_T^i)``. This is the Spectral Theorem; its
statement, its proof, its dual version in frequency space and the link with FinancialFFT.jl are in
[The Spectral Theorem](spectral.md).
