# The Spectral Theorem

In Black and Scholes every sensitivity of order ``n`` of a European option, with respect to any input, is a linear
combination of the prices of ``2n + 1`` contracts, with coefficients that do not depend on the payoff. This is the
Spectral Theorem. It is the meeting point of the two halves of this family of packages:

- VibratoMonteCarlo.jl differentiates the **log density** of the price at maturity, through the likelihood ratio
  weights of [Pathwise and likelihood ratio methods](lrm.md), computed by `v_mod`.
- [FinancialFFT.jl](https://github.com/rcalxrc08/FinancialFFT.jl) differentiates the **log characteristic function**,
  through automatic differentiation of `exp` of the characteristic exponent in the Lewis formula [lewis2001](@cite).

Both are the same Faà di Bruno formula, applied to ``\exp`` of a logarithm, on the two sides of the Fourier transform.
The Spectral Theorem is its statement in price space. Its **dual**, in frequency space, holds well beyond Black and
Scholes: it turns the likelihood ratio weights into derivatives of the payoff, and links the Monte Carlo methods of this
package to the Fourier methods of FinancialFFT.jl.

!!! info "History"
    The differentiation of the FFT methods and the Vibrato Monte Carlo methods were both part of the thesis
    [scaramuzzino2017](@cite). The Spectral Theorem for Black and Scholes came later, in 2024, by chance.
    Its dual version on the complex plane is what ties the two parts together: this page states it, in the polynomial
    form that reduces to derivatives of the payoff and in general, for every parameter of every model with a
    characteristic exponent. Its implementation in FinancialFFT.jl is the next step.

## Setting

The inputs are the spot ``S_0``, the risk free rate ``r``, the dividend yield ``d``, the maturity ``T`` and the
parameters of the model, for Black and Scholes the volatility ``\sigma``. ``\theta`` denotes any one of them, the input
we differentiate with respect to. As in [Pathwise and likelihood ratio methods](lrm.md), ``X_T = \log(S_T)``, the payoff
is ``\Phi(x) = \phi(e^x)``, ``f(x, \theta)`` is the density of ``X_T``, and up to [The general dual theorem](@ref) the
price is not discounted:

```math
V_0 = \mathbb{E}(\Phi(X_T))
```

The discount only adds known terms in ``r``, see the rate in [The dual theorem](@ref dual-spectral-theorem).

The complete exponential Bell polynomials ``B_n(x_1, \dots, x_n)`` (listed in [Bell polynomials up to order 6](@ref))
are the polynomials of Faà di Bruno's formula for the exponential: for a smooth function ``g(\theta)``,

```math
\partial^n_\theta e^{g(\theta)} = e^{g(\theta)}\, B_n(\partial_\theta g, \partial^2_\theta g, \dots, \partial^n_\theta g)
```

With ``g = \log f`` this gives the likelihood ratio formula of order ``n``, whose weight is
``B_n(\{\partial^{i}_{\theta} \log{f}\}_{i = 1 : n})``: the derivatives ``\partial^{i}_{\theta} \log{f}`` are the
scores. Two families of contracts appear:

```math
G_k = \mathbb{E}(\Phi(X_T) X_T^k), \qquad D_k = \mathbb{E}(\Phi^{(k)}(X_T))
```

the moment weighted contracts ``G_k``, which come from the likelihood ratio weights, and the contracts ``D_k`` on the
derivatives of the payoff. These derivatives are taken as distributions: for a call ``\Phi(x) = (e^x - K)^+`` with
strike ``K``, ``\Phi'(x) = e^x \mathbb{1}_{x > \log K}`` and ``\Phi''`` has a Dirac mass ``K \delta_{\log K}``, so ``D_2``
contains the density of ``X_T`` at the strike.

In Black and Scholes, ``X_T = \log S_0 + (r - d - \frac{\sigma^2}{2}) T + \sigma W_T`` with ``W`` a Brownian motion, so
``X_T`` is Gaussian with mean and standard deviation

```math
\mu = \log S_0 + (r - d - \frac{\sigma^2}{2}) T, \qquad s = \sigma \sqrt{T}
```

## The theorem

The density of ``X_T`` in Black and Scholes is Gaussian, so its logarithm is a second order polynomial in ``x``, with
coefficients ``b_0``, ``b_1``, ``b_2`` that depend on the inputs:

```math
\log{f(x,\theta)}= b_0(\theta) + b_1(\theta)\,x +  b_2(\theta)\,x^2, \qquad b_1(\theta) = \frac{\mu}{s^2}, \qquad b_2(\theta) = -\frac{1}{2 s^2} = -\frac{1}{2 \sigma^2 T}
```

Every score is again a polynomial of degree at most 2 in ``x``:

```math
\partial^{i}_{\theta} \log{f(x,\theta)}= \partial^{i}_{\theta}b_0(\theta) + \partial^{i}_{\theta}b_1(\theta)\,x +  \partial^{i}_{\theta}b_2(\theta)\,x^2
```

Each monomial of ``B_n`` is a product of at most ``n`` scores, so the likelihood ratio weight of order ``n`` is a
polynomial of degree at most ``2 n`` in ``x``, whose coefficients we call ``a_i(\theta, n)``:

```math
B_n(\{\partial^{i}_{\theta} \log{f(x,\theta)}\}_{i = 1 : n}) = \sum_{i=0}^{2 n} a_i(\theta,n)\, x^i
```

The likelihood ratio formula ``\partial^{n}_{\theta} V_0 = \mathbb{E}(\Phi(X_T) B_n(\{\partial^{i}_{\theta} \log{f(X_T,\theta)}\}_{i = 1 : n}))`` then becomes:

!!! note "Spectral Theorem"
    Given a European option with payoff ``\Phi`` on an underlying that follows the Black and Scholes model, a
    sensitivity of order ``n`` with respect to any input ``\theta`` is

    ```math
    \partial^{n}_{\theta} V_0  = \sum_{i=0}^{2 n} a_i(\theta,n)\, G_i = a_0(\theta,n) V_0 + \sum_{i=1}^{2 n} a_i(\theta,n)\, \mathbb{E}(\Phi(X_T) X_T^i)
    ```

    The coefficients ``a_i(\theta,n)`` do not depend on the payoff: the prices of the ``2 n + 1`` contracts
    ``G_0, \dots, G_{2n}`` give every sensitivity up to order ``n``, with respect to every input.

Two remarks make the count precise.

- ``b_2`` depends on ``\sigma`` and ``T`` only. For ``\theta \in \{S_0, r, d\}`` its derivatives vanish, the scores
  are polynomials of degree 1 and the weight has degree ``n``: ``n + 1`` contracts are enough.
- The coefficients are the same for every payoff and every strike. For a whole smile the contracts ``G_k`` are needed at
  every strike, and each ``k`` is one Fourier transform over the strikes: ``2n + 1`` transforms in all.

### Why it holds: exponential families

The proof uses only one property of the Gaussian density: it is an exponential family whose sufficient statistics are
powers of ``x``. The theorem holds unchanged for every density of the form

```math
f(x,\theta) = h(x) \exp\Big(\sum_{j=1}^{m} \eta_j(\theta)\, x^j - A(\theta)\Big)
```

with ``h \ge 0`` a function of ``x`` only, ``\eta_1, \dots, \eta_m`` and ``A`` functions of ``\theta`` only. Indeed
``\partial^i_\theta \log f = \sum_j \partial^i_\theta \eta_j(\theta)\, x^j - \partial^i_\theta A(\theta)`` is a polynomial
of degree at most ``m``, and the weight of order ``n`` a polynomial of degree at most ``m n``: ``\partial^n_\theta V_0``
is a combination of ``G_0, \dots, G_{mn}``. Black and Scholes is the case ``m = 2``. The densities of the Lévy models
(Variance Gamma, Normal Inverse Gaussian, Kou, Merton) are not of this form: their logarithms involve Bessel functions or
series, so in price space the theorem does not extend to them. It does extend in frequency space, see below.

### Hermite form and integration by parts

For the inputs that enter only the mean of ``X_T``, the score polynomials are Hermite polynomials. With
``\delta = \sqrt{T}`` and ``Z = (x - \mu) / (\sigma \delta)``, which is a standard normal when ``x = X_T`` (as in
[Likelihood Ratio Method for Black and Scholes](@ref)), Symbolics.jl generates from the `log_density` of the package the
weights for ``\theta = r``:

```@eval
Main.DocumentationFormulas.hermite_r_weights()
```

These are ``(\delta / \sigma)^n He_n(Z)``, with ``He_n`` the probabilists' Hermite polynomials (``He_0 = 1``,
``He_1(z) = z``, ``He_{n+1}(z) = z He_n(z) - n He_{n-1}(z)``); the documentation build checks it numerically. Gaussian
integration by parts, ``\mathbb{E}(g(Z) He_n(Z)) = \mathbb{E}(g^{(n)}(Z))`` for a standard normal ``Z``, applied to
``g(z) = \Phi(\mu + \sigma \delta z)``, whose ``n``-th derivative is ``(\sigma \delta)^n \Phi^{(n)}(\mu + \sigma \delta z)``,
turns the weight into derivatives of the payoff:

```math
\partial^{n}_{r} V_0 = \frac{\delta^n}{\sigma^n} \mathbb{E}(\Phi(X_T) He_n(Z)) = \frac{\delta^n}{\sigma^n} (\sigma \delta)^n \mathbb{E}(\Phi^{(n)}(X_T)) = T^n D_n
```

The same sensitivity is a combination of the moment weighted contracts ``G_k`` and a single contract ``D_n`` on the
derivative of the payoff. This is the first instance of the duality below.

### Application to a sixth order derivative

Let's assume ``\theta = r``. From the Spectral Theorem, in terms of ``W_T = (X_T - \mu) / \sigma``:
```@eval
Main.DocumentationFormulas.sixth_order()
```
The weight has only even powers of ``W_T``, so we need the 4 contracts ``\mathbb{E}(\Phi(X_T) W_T^k)``,
``k = 0, 2, 4, 6``, one of them the option itself. In Hermite form the weight is ``T^3 He_6(W_T / \sqrt{T}) / \sigma^6``,
and by the previous section the sensitivity is also ``T^6 D_6``.

## The Fourier dictionary

Let ``L`` be the logarithm of the moment generating function of ``X_T``, defined for complex ``w`` where the expectation
is finite, the function that FinancialFFT.jl builds from the characteristic exponent of the model
(`FinancialFFT.characteristic_exponent_i`). On the imaginary axis, ``w = i u`` with ``u`` real, it gives the
characteristic function:

```math
L(w) = \log \mathbb{E}(e^{w X_T}), \qquad \mathbb{E}(e^{i u X_T}) = e^{L(i u)}
```

Let ``\hat\Phi(z) = \int_{\mathbb{R}} e^{i z x}\, \Phi(x)\, dx`` be the transform of the payoff, defined for ``z`` in a
strip of the complex plane (``\operatorname{Im} z > 1`` for a call). Inverting ``\hat\Phi`` on a line
``z = \xi + i \beta``, ``\xi \in \mathbb{R}``, inside the strip, where ``\mathbb{E}(e^{\beta X_T}) < \infty``, and taking
the expectation inside the integral gives the Fourier price [lewis2001](@cite), [eberlein2010](@cite):

```math
V_0 = \frac{1}{2 \pi} \int_{\mathbb{R}} \hat\Phi(z)\, e^{L(-i z)}\, d\xi, \qquad z = \xi + i \beta
```

Write ``w = -i z = \beta - i \xi``. The two families of contracts become two operations on ``e^{L}``:

- **Differentiating the payoff is multiplying by ``w``.** Since ``\Phi^{(k)}(x) = \frac{1}{2\pi}\int (-iz)^k e^{-izx}\hat\Phi(z)\,d\xi``,

  ```math
  D_k = \frac{1}{2 \pi} \int_{\mathbb{R}} \hat\Phi(z)\, w^k\, e^{L(w)}\, d\xi
  ```

- **Weighting by ``X_T^k`` is differentiating in ``w``.** Since ``\mathbb{E}(X_T^k e^{w X_T}) = \partial_w^k \mathbb{E}(e^{w X_T})``,

  ```math
  G_k = \frac{1}{2 \pi} \int_{\mathbb{R}} \hat\Phi(z)\, \partial_w^k e^{L(w)}\, d\xi
  ```

The likelihood ratio weights act as differential operators in ``w``, the derivatives of the payoff as polynomials in
``w``.

## [The dual theorem](@id dual-spectral-theorem)

The inputs enter the Fourier price only through ``L``, and Faà di Bruno applies to ``e^{L}`` exactly as it applies to
``e^{\log f}`` in the likelihood ratio method. Where the derivatives can be taken under the integral (the hypotheses are
in [The general dual theorem](@ref)),

```math
\partial^n_\theta V_0 = \frac{1}{2 \pi} \int_{\mathbb{R}} \hat\Phi(z)\, e^{L(w)}\, B_n(\{\partial^{i}_{\theta} L(w)\}_{i = 1 : n})\, d\xi, \qquad w = -i z
```

We call ``B_n(\{\partial^{i}_{\theta} L\})(w)`` the multiplier of order ``n``, the frequency space score: like the
likelihood ratio weight, it depends on the model, the input and the maturity, never on the payoff. When it is a polynomial
in ``w``, the first item of the dictionary turns it into derivatives of the payoff:

!!! note "Dual Spectral Theorem"
    If every ``\partial^{i}_{\theta} L(w)``, ``i \le n``, is a polynomial of degree at most ``m`` in ``w``, then the
    multiplier is a polynomial of degree at most ``m n``, with coefficients ``c_k(\theta, n)``, and

    ```math
    B_n(\{\partial^{i}_{\theta} L(w)\}_{i = 1 : n}) = \sum_{k=0}^{m n} c_k(\theta,n)\, w^k, \qquad \partial^{n}_{\theta} V_0 = \sum_{k=0}^{m n} c_k(\theta,n)\, D_k = \mathbb{E}((P_{\theta,n}(\partial_x)\Phi)(X_T))
    ```

    with ``P_{\theta,n}(w) = \sum_k c_k(\theta,n)\, w^k``: the sensitivity is the price of the payoff transformed by a
    differential operator with constant coefficients, which do not depend on the payoff.

The primal theorem needs a polynomial log density, the dual one only polynomial derivatives of ``L``, which is a much
weaker condition. For ``S_0``, ``r`` and ``d`` it holds for **every** Lévy model. There the log price is
``X_T = \log S_0 + (r - d - \psi(1)) T + Y_T``, with ``Y`` a Lévy process with ``Y_0 = 0`` and
``\psi(w) = \log \mathbb{E}(e^{w Y_1})`` its exponent (`characteristic_exponent_i`); ``\psi(1)`` is the martingale
correction, which makes ``e^{-(r - d) t} S_t`` a martingale. Hence

```math
L(w) = w \left(\log S_0 + (r - d) T - T \psi(1)\right) + T \psi(w)
```

and ``S_0``, ``r``, ``d`` appear only in the coefficient of ``w``. For Black and Scholes, ``\psi(w) = \sigma^2 w^2 / 2``,
Symbolics.jl gives:

```@eval
Main.DocumentationFormulas.dual_symbols()
```

and the multipliers of order 2 and 3 for ``S_0`` and ``\sigma``, written ``B_n^{S_0}`` and ``B_n^{\sigma}``:

```@eval
Main.DocumentationFormulas.dual_multipliers()
```

They read as operators on the payoff:

- **Spot.** ``\partial^i_{S_0} L = (-1)^{i-1} (i-1)!\, w / S_0^i``, whose Bell polynomial is the falling factorial
  ``w (w - 1) \cdots (w - n + 1) / S_0^n``: in every Lévy model
  ```math
  \partial^n_{S_0} V_0 = S_0^{-n}\, \mathbb{E}((\partial_x (\partial_x - 1) \cdots (\partial_x - n + 1) \Phi)(X_T))
  ```
  For a call at ``n = 2``, ``\partial_x(\partial_x - 1)\Phi = K \delta_{\log K}``, so the gamma is
  ``\Gamma = K f(\log K) / S_0^2``, with ``f`` the density of ``X_T``.
- **Rate and dividend.** ``\partial_r L = T w``, ``\partial_d L = -T w``, with no higher derivatives:
  ``\partial^n_r V_0 = T^n D_n`` and ``\partial^n_d V_0 = (-T)^n D_n`` in every Lévy model. With the discount
  ``e^{-rT}`` in ``V_0``, ``\partial_r L`` becomes ``T (w - 1)`` and the operator ``T^n (\partial_x - 1)^n``: at
  ``n = 1`` a Black and Scholes call gives ``\rho = T K e^{-rT} N(d_2)``, with ``N`` the standard normal distribution
  function and ``d_2 = (\log(S_0 / K) + (r - d - \sigma^2 / 2) T) / (\sigma \sqrt{T})``.
- **Volatility of a Brownian component.** With ``q(w) = w^2 - w``, ``\partial_\sigma L = T \sigma q`` and
  ``\partial^2_\sigma L = T q``, and the higher derivatives vanish. The multiplier is ``B_n`` evaluated at
  ``x_1 = T \sigma q``, ``x_2 = T q``, ``x_3 = \dots = x_n = 0``, a polynomial of degree ``2n`` in ``w``; as an operator,
  ``q`` becomes ``\mathcal{A} = \partial_x^2 - \partial_x``, so ``2n + 1`` contracts ``D_0, \dots, D_{2n}`` are needed.
  It holds for Black and Scholes, Merton and Kou, whose ``\sigma`` is the volatility of a separate Brownian motion, not
  for Variance Gamma and Normal Inverse Gaussian, whose ``\sigma`` is inside the subordinated Brownian motion.
- **Maturity.** ``\partial_T L = w (r - d - \psi(1)) + \psi(w)`` does not depend on ``T``, so
  ``\partial^n_T V_0 = \mathbb{E}((\mathcal{L}^n \Phi)(X_T))``, where ``\mathcal{L}`` is the operator that multiplies the
  transform by ``\partial_T L``, i.e. ``\mathcal{L} e^{w x} = \partial_T L(w)\, e^{w x}``: the infinitesimal generator of
  ``X``. This is the Kolmogorov backward equation. In Black and Scholes
  ``\mathcal{L} = \frac{\sigma^2}{2} \partial_x^2 + (r - d - \frac{\sigma^2}{2}) \partial_x`` is a differential operator; in
  a Lévy model with jumps it is an integro-differential operator, and only the general form below applies.
- **Jump intensity.** For a compound Poisson part with intensity ``\lambda``, jump sizes distributed as ``J`` and
  ``c = \mathbb{E}(e^J) - 1``, ``\psi`` contains ``\lambda (\mathbb{E}(e^{w J}) - 1)``, so
  ``\partial_\lambda L = T (\mathbb{E}(e^{w J}) - 1 - c w)`` and higher derivatives vanish. The factor
  ``\mathbb{E}(e^{w J})`` is the shift ``\Phi \mapsto \mathbb{E}(\Phi(\cdot + J))``, so with
  ``\Delta_J \Phi(x) = \mathbb{E}(\Phi(x + J)) - \Phi(x)``
  ```math
  \partial^n_\lambda V_0 = T^n\, \mathbb{E}(((\Delta_J - c\, \partial_x)^n \Phi)(X_T))
  ```
  a combination of prices with ``0, \dots, n`` extra jumps and of derivatives of the payoff. The multiplier is not a
  polynomial, but the identity is still exact.

For ``S_0``, ``r`` and ``d`` the identities are a direct differentiation of ``\Phi(X_T)``, since ``X_T`` depends on them
only through the shift ``\log S_0 + (r - d) T`` and ``Y_T`` does not. Each identity of the list was checked numerically,
against the Fourier price differentiated with `TaylorSeries.Taylor1` up to the third order, on a smooth payoff with the
right side computed in price space by Gauss–Hermite quadrature: Black and Scholes for every input, Merton for ``S_0``,
``\sigma`` and ``\lambda``, all within ``10^{-11}``.

## Where the two theorems meet

In Black and Scholes both theorems hold: the families ``\{G_k\}`` and ``\{D_k\}`` span the same contracts. The reason is
the second item of the dictionary: ``\partial_w e^{L} = L'(w)\, e^{L}`` and, with the mean ``\mu`` and standard deviation
``s`` of [Setting](@ref), ``L(w) = \mu w + s^2 w^2 / 2``, so ``L'(w) = \mu + s^2 w`` is a polynomial and weighting by
``X_T`` is a polynomial in ``w``. At the first order this is Stein's lemma [stein1981](@cite),
``\mathbb{E}(\Phi(X_T) X_T) = \mu\, \mathbb{E}(\Phi(X_T)) + s^2\, \mathbb{E}(\Phi'(X_T))``; at every order, Symbolics.jl
gives:

```@eval
Main.DocumentationFormulas.stein_bridge()
```

Among all laws of ``X_T``, the Gaussian is the only one where this happens. If ``\partial_w e^{L} = Q(w)\, e^{L}`` with
``Q`` a polynomial, then ``L' = Q`` and ``L`` is a polynomial; by Marcinkiewicz's theorem [marcinkiewicz1939](@cite) a
characteristic function ``e^{L(iu)}`` with ``L`` polynomial has ``L`` of degree at most 2, so ``X_T`` is Gaussian.
Outside Black and Scholes the two theorems separate:

| | primal (price space) | dual (frequency space) |
|:--|:--|:--|
| contracts | ``G_k = \mathbb{E}(\Phi(X_T) X_T^k)`` | ``D_k = \mathbb{E}(\Phi^{(k)}(X_T))`` |
| weight | ``B_n(\{\partial^i_\theta \log f\})``, a polynomial in ``x`` | ``B_n(\{\partial^i_\theta L\})``, a polynomial in ``w`` |
| condition | ``\log f`` polynomial in ``x`` | ``\partial^i_\theta L`` polynomial in ``w`` |
| in the packages | `v_mod` on ``\log f``, VibratoMonteCarlo.jl | AD of ``e^{L}`` in the Fourier price, FinancialFFT.jl |
| holds for | Black and Scholes (exponential families) | ``S_0``, ``r``, ``d`` in every Lévy model; ``\sigma`` of a Brownian component; every parameter through the general multiplier |

## The general dual theorem

When ``\partial^i_\theta L`` is not a polynomial (``\kappa`` of Variance Gamma, the jump parameters, ``T`` in a Lévy
model) the multiplier does not reduce to derivatives of the payoff, but it is still explicit, and the Faà di Bruno
formula holds on the complex plane for every parameter of every model with a characteristic exponent. From here on
``V_0`` is the discounted price ``e^{-rT}\, \mathbb{E}(\Phi(X_T))`` and ``L`` carries the discount:

```math
L(w;\theta) = -r T + w \left(\log S_0 + (r - d) T - T \psi(1)\right) + T \psi(w)
```

!!! note "General Dual Spectral Theorem"
    Let ``\Theta`` be an open interval of values of ``\theta`` and ``\beta`` a real number such that ``\hat\Phi`` exists
    on the line ``\operatorname{Im} z = \beta`` and ``\mathbb{E}(e^{\beta X_T}) < \infty`` for every ``\theta \in \Theta``.
    Assume that there is an integrable function ``g`` of ``\xi`` with
    ``|\hat\Phi(\xi + i \beta)\, \partial^k_\theta e^{L(\beta - i \xi;\theta)}| \le g(\xi)`` for every ``\theta \in \Theta``,
    ``\xi \in \mathbb{R}`` and ``k \le n``. Then ``V_0`` is ``n`` times differentiable on ``\Theta`` and, with
    ``z = \xi + i \beta`` and ``w = -i z = \beta - i \xi``,

    ```math
    \partial^n_\theta V_0 = \frac{1}{\pi} \int_0^\infty \operatorname{Re}\left(\hat\Phi(z)\, e^{L(w;\theta)}\, B_n(\{\partial^{i}_{\theta} L(w;\theta)\}_{i = 1 : n})\right) d\xi
    ```

    The multipliers ``B_k(\{\partial^{i}_{\theta} L\})(w)``, ``k \le n``, depend on the model, the parameter and the
    maturity, not on the payoff.

*Proof.* By [The Fourier dictionary](@ref), ``V_0 = \frac{1}{2\pi} \int_{\mathbb{R}} \hat\Phi(z)\, e^{L(w;\theta)}\, d\xi``
(the discount is the term ``-rT`` of ``L``). The bound ``g`` allows to differentiate ``n`` times under the integral
(dominated convergence, one order at a time), and ``\partial^n_\theta e^{L} = e^{L} B_n(\partial_\theta L, \dots, \partial^n_\theta L)``
is Faà di Bruno, as in [Setting](@ref). ``\Phi`` is real, so ``\hat\Phi(-\bar z) = \overline{\hat\Phi(z)}``, and ``X_T``
is real, so ``e^{L(\bar w)} = \overline{e^{L(w)}}``: the integrand at ``-\xi`` is the conjugate of the one at ``\xi``,
and the integral over ``\mathbb{R}`` is twice the integral of the real part over ``\xi > 0``. ∎

The transforms of the payoffs of the package, with strike ``K``, on the lines used:

| payoff | ``\hat\Phi(z)`` | strip | ``\beta`` |
|:--|:--|:--|:--|
| call ``(e^x - K)^+`` | ``K^{i z + 1} / (i z (i z + 1))`` | ``\operatorname{Im} z > 1`` | 1.5 |
| put ``(K - e^x)^+`` | ``K^{i z + 1} / (i z (i z + 1))`` | ``\operatorname{Im} z < 0`` | -0.5 |
| binary call ``\mathbb{1}_{e^x > K}`` | ``-K^{i z} / (i z)`` | ``\operatorname{Im} z > 0`` | 0.5 |
| binary put ``\mathbb{1}_{e^x < K}`` | ``K^{i z} / (i z)`` | ``\operatorname{Im} z < 0`` | -0.5 |

### The multipliers of the models

``S_0``, ``r`` and ``d`` give the polynomial multipliers of [The dual theorem](@ref dual-spectral-theorem) (with the
discount, ``\partial_r L = T(w - 1)``). ``L`` is linear in ``T``, so ``\partial_T L = -r + w (r - d - \psi(1)) + \psi(w)``
and ``B_n = (\partial_T L)^n``. A parameter ``\theta`` of the exponent ``\psi`` enters ``L`` twice, through the law of
``Y`` and through the martingale correction ``\psi(1)``:

```math
\partial^i_\theta L(w) = T \left(\partial^i_\theta \psi(w) - w\, \partial^i_\theta \psi(1)\right)
```

Two consequences hold for every model. The multiplier vanishes at ``w = 0`` and at ``w = 1``, since ``\psi(0) = 0`` for
every value of ``\theta`` and the two terms cancel at ``w = 1``: ``e^{L(0)} = e^{-rT}`` and
``e^{L(1)} = S_0 e^{-dT}``, the prices of the bond and of the forward, do not depend on the parameters of the model. This
is the frequency space counterpart of the zero mean of the likelihood ratio scores. And the correction term grows like
``|w|``, which drives the integrability below.

The exponents ``\psi`` of FinancialFFT.jl (`characteristic_exponent_i`) for the jump and subordinated models are below.
In all of them ``\sigma`` is the volatility of the Brownian part. In Merton, jumps arrive with intensity ``\lambda`` and
are normal with mean ``\mu_J`` and standard deviation ``\sigma_J``. In Kou, jumps arrive with intensity ``\lambda`` and
are exponential, upwards with probability ``p`` and rate ``\lambda_+``, downwards with rate ``\lambda_-``. Variance Gamma
and Normal Inverse Gaussian are a Brownian motion with drift ``\vartheta`` (the field `θ` of the processes) and volatility
``\sigma``, evaluated at a Gamma, respectively Inverse Gaussian, subordinator with mean ``t`` and variance ``\kappa t``.

```math
\begin{aligned}
\psi_{\text{Merton}}(w) &= \frac{\sigma^2 w^2}{2} + \lambda \left(e^{\mu_J w + \sigma_J^2 w^2 / 2} - 1\right) \\
\psi_{\text{Kou}}(w) &= \frac{\sigma^2 w^2}{2} + \lambda w \left(\frac{p}{\lambda_+ - w} - \frac{1 - p}{\lambda_- + w}\right) \\
\psi_{\text{VG}}(w) &= -\frac{1}{\kappa} \log\left(1 - \vartheta \kappa w - \frac{\sigma^2 \kappa w^2}{2}\right) \\
\psi_{\text{NIG}}(w) &= \frac{1}{\kappa} \left(1 - \sqrt{1 - 2 \vartheta \kappa w - \sigma^2 \kappa w^2}\right)
\end{aligned}
```

Symbolics.jl gives their first derivatives in each parameter, the first order multipliers before the correction term
(the documentation build checks each ``\psi`` against `characteristic_exponent_i`):

```@eval
Main.DocumentationFormulas.frequency_scores()
```

Higher orders follow from the same differentiation. The exponents of Merton and Kou are linear in ``\lambda``, so
``B_n = (\partial_\lambda L)^n``, which gives the shift identity of [The dual theorem](@ref dual-spectral-theorem). Mixed
derivatives in two parameters ``\theta_1`` and ``\theta_2`` follow from the multivariate Faà di Bruno formula, for example
``\partial_{\theta_1} \partial_{\theta_2} e^{L} = e^{L} (\partial_{\theta_1} \partial_{\theta_2} L + \partial_{\theta_1} L\, \partial_{\theta_2} L)``,
which is what a `TaylorN` or a nested ForwardDiff dual computes.

### Algorithm and properties

Take ``N_\xi`` nodes ``\xi_j > 0`` of a quadrature of the half line, with weights ``\omega_j``, and the points
``z_j = \xi_j + i \beta``, ``w_j = -i z_j``:

1. evaluate ``L(w_j;\theta)`` with ``\theta`` a `Taylor1` of order ``n``: its coefficients give
   ``\partial^i_\theta L(w_j)``, ``i \le n``;
2. compute the multipliers ``m_{k,j} = e^{L(w_j)} B_k(\{\partial^{i}_{\theta} L(w_j)\})``, ``k \le n``, a matrix of
   nodes × orders (equivalently, ``k!`` times the coefficients of ``\exp`` of the `Taylor1`);
3. for each payoff, ``\partial^k_\theta V_0 \approx \frac{1}{\pi} \sum_j \omega_j \operatorname{Re}(\hat\Phi(z_j)\, m_{k,j})``.

Steps 1 and 2 cost ``O(n^2)`` per node and do not depend on the payoffs; each further payoff, strike or order costs one
weighted sum. A smile with ``N_K`` strikes and every order up to ``n`` is the product of the ``N_K \times N_\xi`` matrix of
weighted payoff transforms ``\omega_j \hat\Phi(z_j)`` with the ``N_\xi \times (n + 1)`` matrix of multipliers. For one
option at one order the method is not cheaper than automatic differentiation through `LewisMethod`, which computes the
same formula implicitly, since ``\exp`` of a truncated Taylor series is Faà di Bruno. What the split adds is the reuse of
the multipliers, the explicit frequency space score (the counterpart of `v_mod`, which does the same for
``\exp(\log f)``), and the reduction to derivatives of the payoff where the multiplier is a polynomial.

The theorem was checked numerically in Julia, with ``S_0 = 100``, ``K = 105``, ``T = 1``, ``r = 0.02``, ``d = 0.01``,
for calls, puts, binary calls and binary puts, on Merton (``S_0``, ``\sigma``, ``\lambda``, ``\mu_J``, ``\sigma_J``), Kou
(``S_0``, ``\sigma``, ``\lambda``, ``p``, ``\lambda_+``, ``\lambda_-``), Variance Gamma and Normal Inverse Gaussian
(``S_0``, ``\sigma``, ``\vartheta``, ``\kappa``, with ``\kappa = 0.2``). Up to the third order it matches `Taylor1` through
`LewisMethod` within ``10^{-9}``, and the first two orders match finite differences of the `LewisMethod` price within
``2 \cdot 10^{-4}`` (the error of the differences). The Black and Scholes call reproduces `blsprice`, its delta and gamma
within ``10^{-14}``.

Differentiating ``\exp`` of ``L`` under the integral is the classical way to get Fourier Greeks: [eberlein2010](@cite)
for the first two orders in ``S_0``, [deolivera2014](@cite) up to the third order, jump parameters included. Its Monte Carlo
counterpart is the likelihood ratio method with a density obtained by inversion [glasserman2010](@cite), which is what
`Vibrato` and `PureLRM` do on Lévy models through the FFT density of FinancialFFT.jl. The form of every order through the
Bell polynomials of ``L``, the payoff independent multipliers and their duality with the Spectral Theorem are the content
of this page; an implementation in FinancialFFT.jl, on the grid of `LewisMethod` and for smiles, is the next step.

### Integrability, and the same limit in Monte Carlo

The integrand behaves like ``|\hat\Phi(z)|\, |e^{L}|\, |B_n|``. With a Brownian component ``|e^{L}|`` decays like a
Gaussian in ``\xi`` and every order converges. In Variance Gamma it decays only like ``|\xi|^{-2T/\kappa}``, while the
multiplier of ``S_0`` grows like ``|\xi|^n`` and ``\hat\Phi`` decays like ``|\xi|^{-2}`` for a call and ``|\xi|^{-1}`` for
a binary. The integral of order ``n`` in ``S_0`` converges absolutely when

```math
n < 1 + \frac{2 T}{\kappa} \quad \text{(call)}, \qquad n < \frac{2 T}{\kappa} \quad \text{(binary)}
```

a sufficient condition, consistent with the failures of Fourier Greeks in Variance Gamma reported by [behrens2023](@cite).
It is the limit met by `VibratoSaltando` on Variance Gamma (see [Methods](methods.md)): with ``a = T / \kappa``, the
``n``-th ``S_0`` derivative of a call at the strike ``K = S_0 e^{(r - d - \psi(1)) T}`` exists only for
``a > (n - 1) / 2``, the same inequality, and one order earlier for a binary. In both cases it is the smoothness of the
true price, seen from the two sides of the transform: the Gamma subordinator ``G_T`` has a density that behaves like
``g^{a - 1}`` at 0, so the density of ``X_T`` is singular at that strike for small ``a``, and its characteristic function
decays slowly.

## The theorem inside the estimators

- `PureLRM` on Black and Scholes computes the primal theorem: its weight of order ``n`` is the polynomial
  ``\sum_i a_i(\theta, n) X_T^i``, evaluated by `v_mod` on the samples, so its derivatives are Monte Carlo estimates of
  ``\sum_i a_i(\theta,n) G_i``.
- `Vibrato` applies it conditionally. Given the path up to ``T - h``, with ``h`` the last time step, the last step is
  Gaussian with a conditional mean ``\mu_h`` and standard deviation ``s_h``, so the derivatives of the inner expectation
  in ``(\mu_h, s_h)`` are combinations of the conditional contracts ``\mathbb{E}(\Phi(X_T) X_T^i \mid \mathcal{F}_{T-h})``
  with deterministic coefficients, ``\mathcal{F}_{T-h}`` being the information up to ``T - h``; the dependence of
  ``\mu_h`` and ``s_h`` on the inputs is differentiated along the simulated path. `VibratoSaltando` does the same given
  the last jump time, or given the subordinator ``G_T`` for Variance Gamma and Normal Inverse Gaussian.
- For ``S_0``, ``r`` and ``d`` the dual theorem turns these conditional weights into derivatives of the payoff, as in
  [Hermite form and integration by parts](@ref).
