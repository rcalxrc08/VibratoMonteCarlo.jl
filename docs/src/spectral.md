# The Spectral Theorem

In Black and Scholes every sensitivity of order ``n`` of a European option, with respect to any parameter, is a linear
combination of the prices of ``2n + 1`` contracts, with coefficients that do not depend on the payoff. This is the
Spectral Theorem. It is the meeting point of the two halves of this family of packages:

- VibratoMonteCarlo.jl differentiates the **log density** of the last step: the likelihood ratio weight
  ``B_n(\{\partial^i_\theta \log f\})`` of [Pathwise and likelihood ratio methods](lrm.md), computed by `v_mod`.
- [FinancialFFT.jl](https://github.com/rcalxrc08/FinancialFFT.jl) differentiates the **log characteristic function**:
  automatic differentiation through `exp` of the characteristic exponent, in the Lewis formula [lewis2001](@cite).

Both are the same Faà di Bruno formula, applied to ``\exp`` of a logarithm, on the two sides of the Fourier transform.
The Spectral Theorem is its statement in price space; its **dual**, in frequency space, holds well beyond Black and
Scholes, and links the likelihood ratio method to the pathwise method and the Monte Carlo methods to the Fourier ones.

!!! info "History"
    The differentiation of the FFT methods and the Vibrato Monte Carlo methods were both part of the thesis
    [scaramuzzino2017](@cite). The Spectral Theorem for Black and Scholes came later, in 2024, by chance.
    Its dual version on the complex plane is what ties the two parts together: this page
    states it, in the polynomial form that reduces to payoff derivatives and in general, for every parameter of every
    model with a characteristic exponent. Its implementation in FinancialFFT.jl is the next step.

## Setting

As in [Pathwise and likelihood ratio methods](lrm.md), ``X_T = \log(S_T)``, the payoff is ``\Phi(x) = \phi(e^x)``, the
density of ``X_T`` is ``f(x, \theta)`` and the price is not discounted:

```math
V_0 = \mathbb{E}(\Phi(X_T))
```

The discount only adds known terms in ``r``, see [The dual theorem](@ref dual-spectral-theorem). Two families of
contracts appear:

```math
G_k = \mathbb{E}(\Phi(X_T) X_T^k), \qquad D_k = \mathbb{E}(\Phi^{(k)}(X_T))
```

the moment weighted contracts ``G_k``, which come from the likelihood ratio method, and the payoff derivative contracts
``D_k``, which come from the pathwise method. The derivatives of ``\Phi`` are taken as distributions: for a call
``\Phi(x) = (e^x - K)^+``, ``\Phi'(x) = e^x \mathbb{1}_{x > \log K}`` and ``\Phi''`` has a Dirac mass ``K \delta_{\log K}``,
so ``D_2`` contains the density of ``X_T`` at the strike. ``B_n`` are the complete exponential Bell polynomials of
[Bell polynomials up to order 6](@ref).

## The theorem

The density of ``X_T`` in Black and Scholes is Gaussian, so its logarithm is a second order polynomial in ``x``, with
coefficients that depend on the parameters:

```math
\log{f(x,\theta)}= b_0(\theta) + b_1(\theta)\,x +  b_2(\theta)\,x^2, \qquad b_2(\theta) = -\frac{1}{2 \sigma^2 T}
```

Every derivative of the log density with respect to ``\theta`` is again a polynomial of degree at most 2 in ``x``:

```math
\partial^{i}_{\theta} \log{f(x,\theta)}= \partial^{i}_{\theta}b_0(\theta) + \partial^{i}_{\theta}b_1(\theta)\,x +  \partial^{i}_{\theta}b_2(\theta)\,x^2
```

``B_n`` has degree ``n`` and is evaluated at polynomials of degree 2, so the likelihood ratio weight of order ``n`` is a
polynomial of degree ``2 n`` in ``x``:

```math
B_n(\{\partial^{i}_{\theta} \log{f(x,\theta)}\}_{i = 1 : n}) = \sum_{i=0}^{2 n} a_i(\theta,n)\, x^i
```

and the likelihood ratio formula ``\partial^{n}_{\theta} V_0 = \mathbb{E}(\Phi(X_T) B_n(\{\partial^{i}_{\theta} \log{f(X_T,\theta)}\}_{i = 1 : n}))`` becomes:

!!! note "Spectral Theorem"
    Given a European option with payoff ``\Phi`` on an underlying that follows the Black and Scholes model, a
    sensitivity of order ``n`` with respect to any parameter ``\theta`` is

    ```math
    \partial^{n}_{\theta} V_0  = \sum_{i=0}^{2 n} a_i(\theta,n)\, G_i = a_0(\theta,n) V_0 + \sum_{i=1}^{2 n} a_i(\theta,n)\, \mathbb{E}(\Phi(X_T) X_T^i)
    ```

    The coefficients ``a_i(\theta,n)`` do not depend on the payoff: the prices of the ``2 n + 1`` contracts
    ``G_0, \dots, G_{2n}`` give every sensitivity up to order ``n``, with respect to every parameter.

Two remarks make the count precise.

- ``b_2`` depends on ``\sigma`` and ``T`` only. For ``\theta \in \{S_0, r, d\}`` its derivatives vanish, the scores
  are polynomials of degree 1 and the weight has degree ``n``: ``n + 1`` contracts are enough.
- The coefficients are the same for every payoff and every strike. The contracts ``G_k`` of a whole smile are ``2n + 1``
  Fourier integrals, one per ``k``, so that a single set of prices gives the sensitivities of every option of the smile.

### Why it holds: exponential families

The proof uses only one property of the Gaussian density: it is an exponential family whose sufficient statistics are
powers of ``x``. The theorem holds unchanged for every density of the form

```math
f(x,\theta) = h(x) \exp\Big(\sum_{j=1}^{p} \eta_j(\theta)\, x^j - A(\theta)\Big)
```

since ``\partial^i_\theta \log f = \sum_j \partial^i_\theta \eta_j(\theta)\, x^j - \partial^i_\theta A(\theta)``
is a polynomial of degree ``p``, and the weight of order ``n`` is a polynomial of degree ``p n``: ``\partial^n_\theta V_0``
is a combination of ``G_0, \dots, G_{pn}``. Black and Scholes is the case ``p = 2``. The densities of the Lévy models
(Variance Gamma, Normal Inverse Gaussian, Kou, Merton) are not of this form: their logarithms involve Bessel functions or
series, so in price space the theorem does not extend to them. It does extend in frequency space, see below.

### Hermite form and integration by parts

For the parameters that enter only the mean of ``X_T``, the score polynomials are Hermite polynomials. With
``\delta = \sqrt{T}`` and ``Z`` the standard normal of [Likelihood Ratio Method for Black and Scholes](@ref), Symbolics.jl
generates from the `log_density` of the package the weights for ``\theta = r``:

```@eval
Main.DocumentationFormulas.hermite_r_weights()
```

These are ``(\delta / \sigma)^n He_n(Z)``, with ``He_n`` the probabilists' Hermite polynomials (the generator checks it).
Gaussian integration by parts, ``\mathbb{E}(g(Z) He_n(Z)) = \mathbb{E}(g^{(n)}(Z))``, applied to
``g(z) = \Phi(\mu + \sigma \delta z)``, turns the weight into derivatives of the payoff:

```math
\partial^{n}_{r} V_0 = \frac{\delta^n}{\sigma^n} \mathbb{E}(\Phi(X_T) He_n(Z)) = \frac{\delta^n}{\sigma^n} (\sigma \delta)^n \mathbb{E}(\Phi^{(n)}(X_T)) = T^n D_n
```

The same sensitivity is a combination of the ``G_k`` (likelihood ratio) and a single ``D_n`` (pathwise). This is the
first instance of the duality below.

### Application to a sixth order derivative

Let's assume ``\theta = r``. From the Spectral Theorem:
```@eval
Main.DocumentationFormulas.sixth_order()
```
Hence we need to price 4 contracts to get the sixth order sensitivity, one of them the option itself. In Hermite form the
weight is ``T^3 He_6(W_T / \sqrt{T}) / \sigma^6``, and by the previous section the sensitivity is also ``T^6 D_6``.

## The Fourier dictionary

Let ``L`` be the logarithm of the moment generating function of ``X_T``, the function that FinancialFFT.jl builds from the
characteristic exponent of the model (`FinancialFFT.characteristic_exponent_i`):

```math
L(w) = \log \mathbb{E}(e^{w X_T}), \qquad \mathbb{E}(e^{i u X_T}) = e^{L(i u)}
```

and ``\hat\Phi(z) = \int_{\mathbb{R}} e^{i z x}\, \Phi(x)\, dx`` the transform of the payoff, defined on a strip of the
complex plane (``\operatorname{Im} z > 1`` for a call). Inverting ``\hat\Phi`` on a line ``\operatorname{Im} z = \beta``
inside the strip, where ``\mathbb{E}(e^{\beta X_T}) < \infty``, and taking the expectation inside the integral gives the
Fourier price [lewis2001](@cite), [eberlein2010](@cite):

```math
V_0 = \frac{1}{2 \pi} \int_{\operatorname{Im} z = \beta} \hat\Phi(z)\, e^{L(-i z)}\, dz
```

Write ``w = -i z``, on the line ``\operatorname{Re} w = \beta``. The two families of contracts become two operations on
``e^{L}``:

- **Differentiating the payoff is multiplying by ``w``.** Since ``\Phi^{(k)}(x) = \frac{1}{2\pi}\int (-iz)^k e^{-izx}\hat\Phi(z)\,dz``,

  ```math
  D_k = \frac{1}{2 \pi} \int \hat\Phi(z)\, w^k\, e^{L(w)}\, dz
  ```

- **Weighting by ``X_T^k`` is differentiating in ``w``.** Since ``\mathbb{E}(X_T^k e^{w X_T}) = \partial_w^k \mathbb{E}(e^{w X_T})``,

  ```math
  G_k = \frac{1}{2 \pi} \int \hat\Phi(z)\, \partial_w^k e^{L(w)}\, dz
  ```

The likelihood ratio method works with differential operators in ``w``, the pathwise method with polynomials in ``w``.

## [The dual theorem](@id dual-spectral-theorem)

The parameters enter the Fourier price only through ``L``, and Faà di Bruno applies to ``e^{L}`` exactly as it applies to
``e^{\log f}`` in the likelihood ratio method:

```math
\partial^n_\theta V_0 = \frac{1}{2 \pi} \int_{\operatorname{Im} z = \beta} \hat\Phi(z)\, e^{L(w)}\, B_n(\{\partial^{i}_{\theta} L(w)\}_{i = 1 : n})\, dz, \qquad w = -i z
```

The multiplier ``B_n(\{\partial^{i}_{\theta} L\})(w)`` is the frequency space score: like the likelihood ratio weight, it
depends on the model, the parameter and the maturity, never on the payoff. When it is a polynomial in ``w``, the first
item of the dictionary turns it into payoff derivatives:

!!! note "Dual Spectral Theorem"
    If every ``\partial^{i}_{\theta} L(w)``, ``i \le n``, is a polynomial of degree at most ``p`` in ``w``, then

    ```math
    B_n(\{\partial^{i}_{\theta} L(w)\}_{i = 1 : n}) = \sum_{k=0}^{p n} c_k(\theta,n)\, w^k, \qquad \partial^{n}_{\theta} V_0 = \sum_{k=0}^{p n} c_k(\theta,n)\, D_k = \mathbb{E}((P_{\theta,n}(\partial_x)\Phi)(X_T))
    ```

    with ``P_{\theta,n}(w) = \sum_k c_k(\theta,n)\, w^k``: the sensitivity is the price of the payoff transformed by a
    differential operator with constant coefficients, which do not depend on the payoff.

The primal theorem needs a polynomial log density, the dual one only polynomial derivatives of ``L``, which is a much
weaker condition. For the model inputs it holds for **every** Lévy model. With the driftless exponent
``\psi(w) = \log \mathbb{E}(e^{w Y_1})`` (`characteristic_exponent_i`), the risk neutral log price has

```math
L(w) = w \left(\log S_0 + (r - d) T - T \psi(1)\right) + T \psi(w)
```

and ``S_0``, ``r``, ``d`` appear only in the coefficient of ``w``. For Black and Scholes, ``\psi(w) = \sigma^2 w^2 / 2``,
Symbolics.jl gives:

```@eval
Main.DocumentationFormulas.dual_symbols()
```

and the multipliers of order 2 and 3 for ``S_0`` and ``\sigma``:

```@eval
Main.DocumentationFormulas.dual_multipliers()
```

They read as operators on the payoff:

- **Spot.** ``\partial^i_{S_0} L = (-1)^{i-1} (i-1)!\, w / S_0^i``, whose Bell polynomial is the falling factorial
  ``w (w - 1) \cdots (w - n + 1) / S_0^n``: in every Lévy model
  ```math
  \partial^n_{S_0} V_0 = S_0^{-n}\, \mathbb{E}((\partial_x (\partial_x - 1) \cdots (\partial_x - n + 1) \Phi)(X_T))
  ```
  For a call at ``n = 2``, ``\partial_x(\partial_x - 1)\Phi = K \delta_{\log K}``, so ``\Gamma = K f(\log K) / S_0^2``,
  the density of ``X_T`` at the strike.
- **Rate and dividend.** ``\partial_r L = T w``, ``\partial_d L = -T w``, with no higher derivatives:
  ``\partial^n_r V_0 = T^n D_n`` and ``\partial^n_d V_0 = (-T)^n D_n`` in every Lévy model. With the discount
  ``e^{-rT}`` in ``V_0``, ``\partial_r L`` becomes ``T (w - 1)`` and the operator ``T^n (\partial_x - 1)^n``: at
  ``n = 1`` a call gives ``\rho = T K e^{-rT} N(d_2)``.
- **Volatility of a Brownian component.** With ``q(w) = w^2 - w``, ``\partial_\sigma L = T \sigma q`` and
  ``\partial^2_\sigma L = T q``: the operator is ``B_n(T \sigma \mathcal{A}, T \mathcal{A})`` with
  ``\mathcal{A} = \partial_x^2 - \partial_x``, so ``2n + 1`` contracts ``D_0, \dots, D_{2n}``. It holds for Black and
  Scholes, Merton and Kou, whose ``\sigma`` is the volatility of a separate Brownian motion, not for Variance Gamma and
  Normal Inverse Gaussian, whose ``\sigma`` is inside the subordinated Brownian motion.
- **Maturity.** ``\partial_T L = w (r - d - \psi(1)) + \psi(w)`` does not depend on ``T``, so
  ``\partial^n_T V_0 = \mathbb{E}((\mathcal{L}^n \Phi)(X_T))`` with ``\mathcal{L}`` the operator of symbol ``\partial_T L``: the
  infinitesimal generator of ``X``. This is the Kolmogorov backward equation. In Black and Scholes
  ``\mathcal{L} = \frac{\sigma^2}{2} \partial_x^2 + (r - d - \frac{\sigma^2}{2}) \partial_x`` is a differential operator; in
  a Lévy model it is an integro-differential operator, and the theorem holds in its general form below.
- **Jump intensity.** For a compound Poisson part with jumps ``J`` and ``c = \mathbb{E}(e^J) - 1``,
  ``\partial_\lambda L = T (\mathbb{E}(e^{w J}) - 1 - c w)`` and higher derivatives vanish. The factor
  ``\mathbb{E}(e^{w J})`` is the shift ``\Phi \mapsto \mathbb{E}(\Phi(\cdot + J))``, so with
  ``\Delta_J \Phi(x) = \mathbb{E}(\Phi(x + J)) - \Phi(x)``
  ```math
  \partial^n_\lambda V_0 = T^n\, \mathbb{E}(((\Delta_J - c\, \partial_x)^n \Phi)(X_T))
  ```
  a combination of prices with ``0, \dots, n`` extra jumps and their payoff derivatives.

For ``S_0``, ``r`` and ``d`` the dual theorem is the pathwise method itself, since ``X_T = \log S_0 + (r - d) T + Y_T``
with ``Y_T`` independent of them. For ``\sigma`` and ``T`` in Black and Scholes it is the pathwise method followed by Gaussian integration by parts. Each
identity of the list was checked numerically, against the Fourier price differentiated with `TaylorSeries.Taylor1` up to
the third order, on a smooth payoff with the right side computed in price space by Gauss–Hermite quadrature: Black and
Scholes for every parameter, Merton for ``S_0``, ``\sigma`` and ``\lambda``, all within ``10^{-11}``.

## Where the two theorems meet

In Black and Scholes both theorems hold: the bases ``\{G_k\}`` and ``\{D_k\}`` span the same contracts. The reason is
the second item of the dictionary: ``\partial_w e^{L} = L'(w)\, e^{L}``, and with
``\mu = \log S_0 + (r - d - \frac{\sigma^2}{2}) T`` and ``s = \sigma \sqrt{T}`` the derivative ``L'(w) = \mu + s^2 w`` is a
polynomial, so weighting by ``X_T`` is a polynomial in ``w``. At the first order this is Stein's lemma
[stein1981](@cite), ``\mathbb{E}(\Phi(X_T) X_T) = \mu\, \mathbb{E}(\Phi(X_T)) + s^2\, \mathbb{E}(\Phi'(X_T))``;
at every order, Symbolics.jl gives:

```@eval
Main.DocumentationFormulas.stein_bridge()
```

Black and Scholes is the only model where this happens, among all laws of ``X_T``. If ``\partial_w e^{L} = q(w)\, e^{L}`` with ``q`` a polynomial,
then ``L' = q`` and ``L`` is a polynomial; by Marcinkiewicz's theorem [marcinkiewicz1939](@cite) a characteristic function
``e^{L(iu)}`` with ``L`` polynomial has degree at most 2, so ``X_T`` is Gaussian. Outside Black and Scholes the two
theorems separate:

| | primal (price space) | dual (frequency space) |
|:--|:--|:--|
| contracts | ``G_k = \mathbb{E}(\Phi(X_T) X_T^k)`` | ``D_k = \mathbb{E}(\Phi^{(k)}(X_T))`` |
| weight | ``B_n(\{\partial^i_\theta \log f\})``, a polynomial in ``x`` | ``B_n(\{\partial^i_\theta L\})``, a polynomial in ``w`` |
| condition | ``\log f`` polynomial in ``x`` | ``\partial^i_\theta L`` polynomial in ``w`` |
| Monte Carlo reading | likelihood ratio | pathwise |
| in the packages | `v_mod` on ``\log f``, VibratoMonteCarlo.jl | AD of ``e^{L}`` in the Fourier price, FinancialFFT.jl |
| holds for | Black and Scholes (exponential families) | ``S_0``, ``r``, ``d`` in every Lévy model; ``\sigma`` of a Brownian component; all parameters through the general multiplier |

## The general dual theorem

When ``\partial^i_\theta L`` is not a polynomial (``\kappa`` of Variance Gamma, the jump parameters, ``T`` in a Lévy
model) the multiplier does not reduce to payoff derivatives, but it is still explicit, and the Faà di Bruno formula holds
on the complex plane for every parameter of every model with a characteristic exponent. From here on ``V_0`` is the
discounted price ``e^{-rT}\, \mathbb{E}(\Phi(X_T))`` and ``L`` carries the discount:

```math
L(w;\theta) = -r T + w \left(\log S_0 + (r - d) T - T \psi(1)\right) + T \psi(w)
```

!!! note "General Dual Spectral Theorem"
    Let ``\Theta`` be an open interval and ``\beta`` a real number such that ``\hat\Phi`` exists on the line
    ``\operatorname{Im} z = \beta`` and ``\mathbb{E}(e^{\beta X_T}) < \infty`` for every ``\theta \in \Theta``. Assume
    that there is an integrable ``g`` with ``|\hat\Phi(\xi + i \beta)\, \partial^k_\theta e^{L(\beta - i \xi;\theta)}| \le g(\xi)``
    for every ``\theta \in \Theta``, ``\xi \in \mathbb{R}`` and ``k \le n``. Then ``V_0`` is ``n`` times differentiable
    on ``\Theta`` and, with ``z = \xi + i \beta`` and ``w = -i z = \beta - i \xi``,

    ```math
    \partial^n_\theta V_0 = \frac{1}{\pi} \int_0^\infty \operatorname{Re}\left(\hat\Phi(z)\, e^{L(w;\theta)}\, B_n(\{\partial^{i}_{\theta} L(w;\theta)\}_{i = 1 : n})\right) d\xi
    ```

    The multipliers ``B_k(\{\partial^{i}_{\theta} L\})(w)``, ``k \le n``, depend on the model, the parameter and the
    maturity, not on the payoff.

*Proof.* By [The Fourier dictionary](@ref), ``V_0 = \frac{1}{2\pi} \int_{\mathbb{R}} \hat\Phi(z)\, e^{L(w;\theta)}\, d\xi``.
The bound ``g`` allows to differentiate ``n`` times under the integral (dominated convergence, one order at a time), and
``\partial^n_\theta e^{L} = e^{L} B_n(\partial_\theta L, \dots, \partial^n_\theta L)`` is Faà di Bruno for ``\exp``, as in
[Likelihood Ratio Method](@ref). ``\Phi`` is real, so ``\hat\Phi(-\bar z) = \overline{\hat\Phi(z)}``, and ``X_T`` is real,
so ``e^{L(\bar w)} = \overline{e^{L(w)}}``: the integrand at ``-\xi`` is the conjugate of the one at ``\xi``, and the
integral over ``\mathbb{R}`` is twice the integral of the real part over ``\xi > 0``. ∎

The transforms of the payoffs of the package, on the lines used:

| payoff | ``\hat\Phi(z)`` | strip | ``\beta`` |
|:--|:--|:--|:--|
| call ``(e^x - K)^+`` | ``K^{i z + 1} / (i z (i z + 1))`` | ``\operatorname{Im} z > 1`` | 1.5 |
| put ``(K - e^x)^+`` | ``K^{i z + 1} / (i z (i z + 1))`` | ``\operatorname{Im} z < 0`` | -0.5 |
| binary call ``\mathbb{1}_{e^x > K}`` | ``-K^{i z} / (i z)`` | ``\operatorname{Im} z > 0`` | 0.5 |
| binary put ``\mathbb{1}_{e^x < K}`` | ``K^{i z} / (i z)`` | ``\operatorname{Im} z < 0`` | -0.5 |

### The multipliers of the models

``S_0``, ``r`` and ``d`` give the polynomial multipliers of [The dual theorem](@ref dual-spectral-theorem) (with the
discount, ``\partial_r L = T(w - 1)``). ``L`` is linear in ``T``, so ``\partial_T L = -r + w (r - d - \psi(1)) + \psi(w)``
and ``B_n = (\partial_T L)^n``. A parameter ``\theta`` of the exponent enters ``L`` twice, through the law of the
increments and through the martingale correction:

```math
\partial^i_\theta L(w) = T \left(\partial^i_\theta \psi(w) - w\, \partial^i_\theta \psi(1)\right)
```

Two consequences hold for every model. The multiplier vanishes at ``w = 0`` and at ``w = 1``, since ``\psi(0) = 0`` for
every ``\theta`` and the correction cancels at ``w = 1``: the bond (``\Phi = 1``) and the forward (``\Phi = e^x``) do not
depend on the parameters of the model, the frequency space counterpart of the zero mean of the likelihood ratio score.
And the correction term grows like ``|w|``, which drives the integrability below.

The exponents ``\psi`` of FinancialFFT.jl (`characteristic_exponent_i`) are, for Merton, Kou, Variance Gamma and Normal
Inverse Gaussian,

```math
\begin{aligned}
\psi_{\text{Merton}}(w) &= \frac{\sigma^2 w^2}{2} + \lambda \left(e^{\mu_J w + \sigma_J^2 w^2 / 2} - 1\right) \\
\psi_{\text{Kou}}(w) &= \frac{\sigma^2 w^2}{2} + \lambda w \left(\frac{p}{\lambda_+ - w} - \frac{1 - p}{\lambda_- + w}\right) \\
\psi_{\text{VG}}(w) &= -\frac{1}{\kappa} \log\left(1 - \theta \kappa w - \frac{\sigma^2 \kappa w^2}{2}\right) \\
\psi_{\text{NIG}}(w) &= \frac{1}{\kappa} \left(1 - \sqrt{1 - 2 \theta \kappa w - \sigma^2 \kappa w^2}\right)
\end{aligned}
```

and Symbolics.jl gives their first derivatives, the first order frequency scores before the correction term (the
generator checks each ``\psi`` against `characteristic_exponent_i`):

```@eval
Main.DocumentationFormulas.frequency_scores()
```

Higher orders follow from the same differentiation. The exponents of Merton and Kou are linear in ``\lambda``, so
``B_n = (\partial_\lambda L)^n``, which gives the shift identity of [The dual theorem](@ref dual-spectral-theorem). Mixed
derivatives follow from the multivariate Faà di Bruno formula, for example
``\partial_{\theta_1} \partial_{\theta_2} e^{L} = e^{L} (\partial_{\theta_1} \partial_{\theta_2} L + \partial_{\theta_1} L\, \partial_{\theta_2} L)``,
which is what a `TaylorN` or a nested ForwardDiff dual computes.

### Algorithm and properties

On the nodes ``\xi_j`` of a quadrature of the line, with weights ``\omega_j``:

1. evaluate ``L(w_j;\theta)`` with ``\theta`` a `Taylor1` of order ``n``: its coefficients give
   ``\partial^i_\theta L(w_j)``, ``i \le n``;
2. compute the multipliers ``m_{k,j} = e^{L(w_j)} B_k(\{\partial^{i}_{\theta} L(w_j)\})``, ``k \le n``, a matrix of
   nodes × orders (equivalently, ``k!`` times the coefficients of ``\exp`` of the `Taylor1`);
3. for each payoff, ``\partial^k_\theta V_0 \approx \frac{1}{\pi} \sum_j \omega_j \operatorname{Re}(\hat\Phi(z_j)\, m_{k,j})``.

Steps 1 and 2 cost ``O(n^2)`` per node and do not depend on the payoffs; each further payoff, strike or order costs one
weighted sum. A smile with ``N_K`` strikes and every order up to ``n`` is the product of the ``N_K \times N_\xi`` matrix of
payoff transforms with the ``N_\xi \times (n + 1)`` matrix of multipliers. For one option at one order the method is not
cheaper than automatic differentiation through `LewisMethod`, which computes the same formula implicitly, since ``\exp``
of a truncated Taylor series is Faà di Bruno. What the split adds is the reuse of the multipliers, the explicit frequency
score (the frequency space counterpart of `v_mod`, which does the same for ``\exp(\log f)``), and the reduction to payoff
derivatives where the multiplier is a polynomial.

The theorem was checked numerically in Julia, with ``S_0 = 100``, ``K = 105``, ``T = 1``, ``r = 0.02``, ``d = 0.01``,
for calls, puts, binary calls and binary puts, on Merton (``S_0``, ``\sigma``, ``\lambda``, ``\mu_J``, ``\sigma_J``), Kou
(``S_0``, ``\sigma``, ``\lambda``, ``p``, ``\lambda_+``, ``\lambda_-``), Variance Gamma and Normal Inverse Gaussian
(``S_0``, ``\sigma``, ``\theta``, ``\kappa``, with ``\kappa = 0.2``). Up to the third order it matches `Taylor1` through
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
Gaussian and every order converges. In Variance Gamma it decays only like ``|z|^{-2T/\kappa}``, while the multiplier of
``S_0`` grows like ``|z|^n`` and ``\hat\Phi`` decays like ``|z|^{-2}`` for a call and ``|z|^{-1}`` for a binary. The
integral of order ``n`` converges absolutely when

```math
n < 1 + \frac{2 T}{\kappa} \quad \text{(call)}, \qquad n < \frac{2 T}{\kappa} \quad \text{(binary)}
```

a sufficient condition, consistent with the failures of Fourier Greeks in Variance Gamma reported by [behrens2023](@cite).
It is the limit met by `VibratoSaltando` on Variance Gamma (see [Methods](methods.md)): with ``a = T / \kappa``, the
``n``-th ``S_0`` derivative of a call exists for ``a > (n - 1) / 2``, the same inequality, and one order earlier for a
binary. In both cases it is the smoothness of the true price, seen from the two sides of the transform: the slow decay of
the characteristic function is the singularity of the density of ``G_T`` at 0.

## The theorem inside the estimators

- `PureLRM` on Black and Scholes computes the primal theorem: the weight of order ``n`` is the polynomial
  ``\sum_i a_i(\theta, n) X_T^i``, evaluated by `v_mod` on the samples, so its derivatives are Monte Carlo estimates of
  ``\sum_i a_i(\theta,n) G_i``.
- `Vibrato` applies it conditionally. Given the path up to ``T - h``, the last step is Gaussian with mean ``\mu`` and
  standard deviation ``s``, so the derivatives of the inner expectation in ``(\mu, s)`` are combinations of the
  conditional contracts ``\mathbb{E}(\Phi(X_T) X_T^i \mid \mathcal{F}_{T-h})`` with deterministic coefficients; the outer
  dependence of ``\mu`` and ``s`` on the parameters is pathwise. `VibratoSaltando` does the same given the last jump
  time, or given ``G_T`` for Variance Gamma and Normal Inverse Gaussian.
- The dual theorem says which of these pieces are pathwise in disguise: for ``S_0``, ``r`` and ``d`` the conditional
  weight integrates to payoff derivatives, as in [Hermite form and integration by parts](@ref).
