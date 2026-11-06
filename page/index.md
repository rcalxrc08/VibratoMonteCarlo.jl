<!-- =============================
     ABOUT
    ============================== -->
\begin{section}{title="About this Package", name="About"}
\lead{VibratoMonteCarlo.jl is a Julia package containing some useful Financial functions for sensitivity computation.}
VibratoMonteCarlo.jl is built on top of [FinancialMonteCarlo.jl](https://github.com/rcalxrc08/FinancialMonteCarlo.jl) and [FinancialFFT.jl](https://github.com/rcalxrc08/FinancialFFT.jl).
Standard montecarlo methods lacks of differentiability, which makes automatic differentiation useless.
The main aim of this package is to provide a feasible way to compute sensitivities of any order for various types of payoffs using montecarlo methods.

The implementation for Vibrato Montecarlo is heavily based on the Master Thesis of Nicola Scaramuzzino \cite{scaramuzzino2017} (that is me) and the subsequent development.

\end{section}
\begin{section}{title="The basic settings"}
Let's assume we have an underlying stock price varying as a stochastic process called $S_t$ and
let's assume without loss of generality that $S_t>0$ a.e..\\
In that case we can define as $X_t$ the following:

$$X_t=\log(S_t) $$

Let's assume we want to price and compute the sensitivities of a european option with maturity $T$ and payoff $\phi$.

We can express the price of such option as:

$$V_0=\mathbb{E}(\phi(S_T)) $$

or, in function of $X_T$:
$$V_0=\mathbb{E}(\Phi(X_T)) $$

Where:

$$\Phi(x)=\phi(e^x) $$
\end{section}
\begin{section}{title="Pathwise Method"}
The pathwise method and the likelihood ratio method are the two classic direct estimators of sensitivities \cite{broadie1996}, see also \cite{glasserman2004}.

Since:
$$V_0=\mathbb{E}(\Phi(X_T))$$

$$\partial_\theta V_0=\partial_\theta \mathbb{E}(\Phi(X_T)) = \mathbb{E}(\dfrac{\partial \Phi(X_T)}{\partial X_T}\dfrac{\partial X_T}{\partial \theta})$$

Where the following is called tangent process:

$$Y_t=\dfrac{\partial X_t}{\partial \theta}$$

In case of a European Option with maturity $T$ and strike price $K$:
$$\Phi(x) = max(e^x-K,0)$$

and:

$$\dfrac{\partial \Phi(x)}{\partial x} = e^x \mathbb{1}_{ e^{x} > K } (x) $$

Hence:

$$\partial_\theta V_0=\mathbb{E}(\dfrac{\partial \Phi(X_T)}{\partial X_T}\dfrac{\partial X_T}{\partial \theta})$$

Unfortunately this method does not provide usefull results for binary options or for n-th order sensitivities, indeed in case of binary options we have:

$$\Phi(x)=\mathbb{1}(x)$$

Hence:

$$\partial_x \Phi(x)=0$$

\end{section}
\begin{section}{title="Likelihood Ratio Method"}
Starting from the price of the basic settings, if we develop the integral:

$$V_0=\mathbb{E}(\Phi(X_T))=\int_{\mathbb{R}}  \Phi(x) f(x) dx $$

where the following is the density of the log yield:

$$f(x)=f_{X_T}(x)$$

If we differentiate for a generic parameter $\theta$:

$$\partial_\theta V_0=\partial_\theta \int_{\mathbb{R}}  \Phi(x) f(x,\theta) dx = \int_{\mathbb{R}}  \Phi(x) \partial_\theta f(x,\theta) dx$$

or:

$$\partial_\theta V_0 = \int_{\mathbb{R}}  \Phi(x) \partial_\theta(\log(f(x,\theta))) f(x,\theta) dx$$

or:

$$\partial_\theta V_0 = \mathbb{E}(\Phi(X_T) \partial_\theta \log(f(X_T,\theta)))$$

In case of n-th order sensitivities:

We observe that for a generic positive differentiable function f we can have:

$$ g(x,\theta) =  \log{f(x,\theta)}$$

or:

$$ f(x,\theta) =  e^{g(x,\theta)}$$

Then, according to [Faa di Bruno's formula](https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula):

$$\partial^{n}_{\theta} e^{g(x,\theta)} = e^{g(x,\theta)} B_n(\{\partial^{i}_{\theta} g(x,\theta)\}_{i = 1 : n})  $$

Where $B_n$ is the complete exponential [Bell polynomial](https://en.wikipedia.org/wiki/Bell_polynomials) of n-th order.
\begin{expand}{Bell polynomials up to order 6}
```julia:bell_polynomials
#hideall
using Symbolics,Latexify,VibratoMonteCarlo
# From the recurrence B_{n+1} = sum_{k=0}^n binomial(n, k) B_{n-k} x_{k+1}.
function generate_bell_polynomials(max_order)
    x_vars = Symbolics.variables(:x, 1:max_order)
    bell = Num[1]
    for n in 0:(max_order-1)
        push!(bell, expand(sum(binomial(n, k) * bell[n-k+1] * x_vars[k+1] for k in 0:n)))
    end
    return bell, x_vars
end
# Terms with more factors first, as they are usually written.
function latex_bell(poly, x_vars)
    coeffs = Symbolics.polynomial_coeffs(poly, x_vars)[1]
    monomials = [(Int(Symbolics.value(c)), [Symbolics.degree(Num(m), x) for x in x_vars]) for (m, c) in coeffs]
    sort!(monomials, by = ((c, e),) -> (-sum(e), findlast(>(0), e)))
    term(c, e) = (c == 1 ? "" : "$c ") * join(("x_$i" * (p > 1 ? "^$p" : "") for (i, p) in enumerate(e) if p > 0), " ")
    return join((term(c, e) for (c, e) in monomials), " + ")
end
bell, x_vars = generate_bell_polynomials(6)
println("\$\$B_0 = 1\$\$")
for n in 1:6
    println("\$\$B_$n(", join(("x_$i" for i in 1:n), ","), ") = ", latex_bell(bell[n+1], x_vars[1:n]), "\$\$")
end
```
\textoutput{bell_polynomials}
\end{expand}

Hence:

$$\partial^{n}_{\theta} f(x,\theta) = f(x,\theta) B_n(\{\partial^{i}_{\theta} \log{f(x,\theta)}\}_{i = 1 : n})  $$

By applying this to the original problem:

$$\partial^{n}_{\theta} V_0=\partial^{n}_{\theta} \int_{\mathbb{R}}  \Phi(x) f(x,\theta) dx = \int_{\mathbb{R}}  \Phi(x) \partial^{n}_{\theta} f(x,\theta) dx$$

By applying the previous result:

$$\partial^{n}_{\theta} V_0= \int_{\mathbb{R}}  \Phi(x) f(x,\theta) B_n(\{\partial^{i}_{\theta} \log{f(x,\theta)}\}_{i = 1 : n}) dx$$

or, in terms of expectation:


$$\partial^{n}_{\theta} V_0 = \mathbb{E}(\Phi(X_T) B_n(\{\partial^{i}_{\theta} \log{f(X_T,\theta)}\}_{i = 1 : n}))$$




\begin{expand}{Likelihood Ratio Method for Black and Scholes}
In case of trivial models the density is known in analytic form, otherwise one need to compute it numerically from the characteristic function.
In case of the Black and Scholes Model \cite{blackscholes1973}, and a European Option we have the following:

$$X_t=\log(S_0)+(r-d-\dfrac{\sigma^2}{2}) t + \sigma W_t$$

$$\log(f(x,S_0,r,d,\sigma,T)) = -\dfrac{1}{2}((\dfrac{x-(\log(S_0)+(r-d-\dfrac{\sigma^2}{2}) T)}{\sigma \sqrt{T}})^2 + \log(2 \pi)) - \log(\sigma\sqrt{T})$$

```julia:lrm_scores
#hideall
# The score is evaluated at x = X_T before simplifying: what is left is a small polynomial in W_T,
# while simplifying the derivative in x first takes hours with SymbolicUtils 4.
function generate_faa_di_bruno(z_better,d_tgt_var,x,X_T)
    init_val=exp(z_better)
    score=expand_derivatives(d_tgt_var(init_val))/init_val
    return simplify_fractions(expand(substitute(score,Dict(x=>X_T))))
end

function generate_faa_di_bruno_exp(z_better,tgt_var,order,x,X_T)
    d_tgt_var=Differential(tgt_var)^order
    return generate_faa_di_bruno(z_better,d_tgt_var,x,X_T)
end

function generate_weight(mu, s, tgt_var, order, x, Z)
    z = VibratoMonteCarlo.log_density(x, mu, s)
    return generate_faa_di_bruno_exp(z, tgt_var, order, x, mu + s * Z)
end
latex_weight(weight) = latexify(weight, env = :raw, mult_symbol = "")
latex_weight(mu, s, tgt_var, order, x, Z) = latex_weight(generate_weight(mu, s, tgt_var, order, x, Z))

function generate_and_adjust(der_sigma,mu,sigma,S_0,r,d,σ,T)
    first_dict=Dict(mu=>log(S_0)+(r-d-σ^2/2)*T,sigma=>σ*sqrt(T));
    new_der_sigma=substitute(der_sigma,first_dict)
	return new_der_sigma
end

function generate_faa_di_bruno_exp_mix(z_better,tgt_var,order,tgt_var2,order2,x,X_T)
    d_tgt_var=Differential(tgt_var)^order *Differential(tgt_var2)^order2
    return generate_faa_di_bruno(z_better,d_tgt_var,x,X_T)
end

function generate_latex_expression_base(S_0,r,d,σ,T,tgt_var,max_order)
    @variables mu,sigma,x
	z=VibratoMonteCarlo.log_density(x, mu, sigma)
    @variables W_T [latexwrapper = identity]
    X_T=log(S_0)+(r-d-σ^2/2)*T+σ*W_T
    z_better=generate_and_adjust(z,mu,sigma,S_0,r,d,σ,T)
    adjusted_final_der=generate_faa_di_bruno_exp(z_better,tgt_var,max_order,x,X_T)
    # der_num=∂ᵣ^max_order*V
    @variables V ∂ ∂θ [latexwrapper = identity]
    der_num=∂^max_order*V / (∂θ^max_order)
    @variables Φ(X_T)
    ress_=Φ*adjusted_final_der
    @variables 𝔼(ress_)
    res=latexify(der_num ~ 𝔼, env=:equation, mult_symbol="")
	return res.s
end
function generate_latex_expression_base_mix(S_0,r,d,σ,T,tgt_var,max_order,tgt_var2,max_order2)
    @variables mu,sigma,x
	z=VibratoMonteCarlo.log_density(x, mu, sigma)
    @variables W_T [latexwrapper = identity]
    order_max=max_order+max_order2
    X_T=log(S_0)+(r-d-σ^2/2)*T+σ*W_T
    z_better=generate_and_adjust(z,mu,sigma,S_0,r,d,σ,T)
    adjusted_final_der=generate_faa_di_bruno_exp_mix(z_better,tgt_var,max_order,tgt_var2,max_order2,x,X_T)
    @variables ∂θ₁ [latexwrapper = identity] V ∂ ∂θ₂ [latexwrapper = identity]
    # der_num=∂ᵣ^max_order*V
    der_num=∂^order_max*V / (∂θ₁^max_order * ∂θ₂^max_order2)
    @variables Φ(X_T)
    ress_=Φ*adjusted_final_der
    @variables 𝔼(ress_)
    res=latexify(der_num ~ 𝔼, env=:equation, mult_symbol="")
	return res.s
end
@variables x S_0 [latexwrapper = identity] r d σ T δ Z
mu_lrm = log(S_0) + (r - d - σ^2 / 2) * δ^2
s_lrm = σ * δ
for (name, tgt_var) in (("S_0", S_0), ("r", r), ("d", d), ("\\sigma", σ))
    println("\$\$\\partial_{", name, "} \\log f = ", latex_weight(mu_lrm, s_lrm, tgt_var, 1, x, Z), "\$\$")
end
# T only enters through δ = sqrt(T), and dδ/dT = 1 / (2δ).
score_T = simplify_fractions(generate_weight(mu_lrm, s_lrm, δ, 1, x, Z) / (2 * δ))
println("\$\$\\partial_{T} \\log f = ", latex_weight(score_T), "\$\$")
```
From this we can compute the various partial derivatives. Calling $\delta = \sqrt{T}$ and:
$$Z = \frac{x - \log(S_0) - (r-d-\frac{\sigma^2}{2}) T}{\sigma \delta}$$
which is a standard normal when $x = X_T$, [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) generates from the `log_density` of the package:
\textoutput{lrm_scores}

Hence, a delta for a generic european option can be computed as the following expectation:

```julia:lrm_delta
#hideall
println("\$\$\\Delta = \\partial_{S_0} V_0 = \\mathbb{E}\\left(\\Phi(X_T) ", latex_weight(mu_lrm, s_lrm, S_0, 1, x, Z), "\\right)\$\$")
```
\textoutput{lrm_delta}

or, in terms of $W_T$:

```julia:lrm_delta_wt
#hideall
println(replace(generate_latex_expression_base(S_0, r, d, σ, T, S_0, 1), "V \\partial" => "\\partial V"))
```
\textoutput{lrm_delta_wt}

Let's assume now that $\Phi$ is a forward i.e.:
$$\Phi(x)=e^x$$

Since $\Delta = e^{(r-d)T}$:
$$\mathbb{E}(e^{-\frac{\sigma^2}{2} T + \sigma W_T}\,W_T) = \sigma T$$
\end{expand}

\begin{expand}{Spectral Theorem for Black and Scholes}
Let's observe something, $\log(f(x,\theta))$ is a second order polymial in $x$:

$$\log{f(X_T,\theta)}= \sum_{i=0}^2 b_i(\theta) X_T^i = b_0(\theta) + b_1(\theta)\,X_T +  b_2(\theta)\,X_T^2$$

Then all of the partial derivatives will be second order polynomials as well:

$$\partial^{i}_{\theta} \log{f(X_T,\theta)}= \partial^{i}_{\theta}b_0(\theta) + \partial^{i}_{\theta}b_1(\theta)\,X_T +  \partial^{i}_{\theta}b_2(\theta)\,X_T^2$$

And the Bell exponential polymial will be a polynomial of degree $2 n$ in terms of $X_T$:

$$B_n(\{\partial^{i}_{\theta} \log{f(X_T,\theta)}\}_{i = 1 : n}) = \sum_{i=0}^{2 n} a_i(\theta,n) X_T^i $$

And finally:

$$\partial^{n}_{\theta} V_0  = a_0(\theta,n) V_0 + \sum_{i=1}^{2 n} a_i(\theta,n) \mathbb{E}(\Phi(X_T) X_T^i) $$

\alert{**Spectral Theorem**: Given a European Option with payoff $\Phi$, an underlying varying according to Black and Scholes model, then we can express a generic $n$ order sensitivity as follows:
$$\partial^{n}_{\theta} V_0  = a_0(\theta,n) V_0 + \sum_{i=1}^{2 n} a_i(\theta,n) \mathbb{E}(\Phi(X_T) X_T^i) $$

Hence a generic $n$ order sensitivity can be expressed as a linear combination of the prices of $2 n + 1$ contracts.
}
To be noticed that the various coefficients $a_i(\theta,n)$ do not depend on the payoff, hence with $2n + 1$ prices we are able to compute all of the partial derivatives with respect to any parameter up to order $n$.

In case $b_2(\theta)$ doesn't depend $\theta$ then the corresponding $n$ order sensitivity can be expressed as a linear combination of the prices of $n + 1$ contracts.
\end{expand}

\begin{expand}{Application for sixth order derivative}
```julia:pure_derivative
#hideall
@variables S_0 [latexwrapper = identity] r d σ T
tgt_var=r
tgt_var2=r
max_order=2
max_order2=4
expr_el=generate_latex_expression_base(S_0,r,d,σ,T,tgt_var,max_order+max_order2)
# expr_el=generate_latex_expression_base_mix(S_0,r,d,σ,T,tgt_var,max_order,tgt_var2,max_order2)
println(expr_el)
```
```julia:mixed_derivative
#hideall
@variables S_0 [latexwrapper = identity] r d σ T
tgt_var=r
tgt_var2=σ
max_order=2
max_order2=1
# expr_el=generate_latex_expression_base(S_0,r,d,σ,T,tgt_var,max_order+max_order2)
expr_el=generate_latex_expression_base_mix(S_0,r,d,σ,T,tgt_var,max_order,tgt_var2,max_order2)
println(expr_el)
```
Let's assume $\theta = r$, from the Spectral Theorem:
\textoutput{pure_derivative}
Hence we need to price 4 contracts to get the sixth order sensitivitity. To be noticed that one of the contract to be priced is the option itself.
\end{expand}

\end{section}
\begin{section}{title="Conditional Expectation"}
The payoff is not smooth, but its conditional expectation is. Let's fix $\tau < T$, let's call $h = T - \tau$ and let's assume that the underlying process is a Levy process. Since:
$$V_0=\mathbb{E}(\Phi(X_T))=\mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))$$
$$\partial^{n}_{\theta} V_0=\partial^{n}_{\theta} \mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))=\mathbb{E}(\partial^{n}_{\theta} \mathbb{E}(\Phi(X_{T})|X_{\tau}))$$
The inner expectation is a smooth function of $X_{\tau}$ and of the parameters, even if $\Phi$ is a kink or a jump. Hence, when it is known in closed form, automatic differentiation can go through it path by path, at any order.

\begin{expand}{Black and Scholes}
The inner expectation is the price of the same option over the remaining time $h$, for a call:
$$C(S,h) = S e^{-d h} N(d_1) - K e^{-r h} N(d_2)$$
$$d_1 = \frac{\log(\frac{S}{K}) + (r-d+\frac{\sigma^2}{2}) h}{\sigma \sqrt{h}}, \qquad d_2 = d_1 - \sigma \sqrt{h}$$
where $N$ is the standard normal distribution function and $\varphi$ its density. Discounting up to $\tau$, the method `Conditional` computes:
$$V_0 = e^{-r \tau}\, \mathbb{E}(C(S_{\tau},h))$$
```julia:conditional_derivatives
#hideall
N(x) = 0.0
φ(x) = 0.0
@register_symbolic N(x)
@register_symbolic φ(x)
@register_derivative N(x) 1 φ(x)
@variables S K h [latexwrapper = identity]
@variables D1(S) D2(S)
@variables d_1 [latexwrapper = identity] d_2 [latexwrapper = identity]
# d_1 and d_2 have the same derivative in S.
function generate_delta(price)
    d_prime = 1 / (S * σ * sqrt(h))
    delta = substitute(expand_derivatives(Differential(S)(price)), Dict(Differential(S)(D1) => d_prime, Differential(S)(D2) => d_prime))
    # The identity S e^{-d h} φ(d_1) = K e^{-r h} φ(d_2).
    return simplify(substitute(delta, Dict(φ(D2) => S * exp(-d * h) * φ(D1) / (K * exp(-r * h)))))
end
latex_delta(delta) = latexify(substitute(delta, Dict(D1 => d_1, D2 => d_2)), env = :raw, mult_symbol = "")
delta_call = generate_delta(S * exp(-d * h) * N(D1) - K * exp(-r * h) * N(D2))
delta_binary = substitute(expand_derivatives(Differential(S)(exp(-r * h) * N(D2))), Dict(Differential(S)(D2) => 1 / (S * σ * sqrt(h))))
println("\$\$\\partial_{S} C_{\\text{call}} = ", latex_delta(delta_call), "\$\$")
println("\$\$\\partial_{S} C_{\\text{binary}} = ", latex_delta(delta_binary), "\$\$")
```
The derivatives of the inner price with respect to $S$, generated with Symbolics, for a call and for a binary call:
\textoutput{conditional_derivatives}
Since $\partial_{S_0} S_{\tau} = \frac{S_{\tau}}{S_0}$, the delta is then:
$$\Delta = e^{-r \tau}\, \mathbb{E}(\frac{S_{\tau}}{S_0} \partial_{S} C(S_{\tau},h))$$
This is the best one can do when it applies, but it needs a closed form for the inner expectation, so it is limited to European and binary options on Gaussian steps. Vibrato removes this limit by computing the same inner expectation by integration, with no closed form.
\end{expand}

\end{section}
\begin{section}{title="Conditional Saltando Expectation"}
Let's assume now that $X_t$ is a finite activity Levy process, like Kou \cite{kou2002} or Merton \cite{merton1976}, with jump intensity $\lambda$ and $\kappa = \mathbb{E}(e^{Y_1} - 1)$, where $Y_i$ are the jump sizes. Instead of a fixed $\tau$, let's take the time of the last jump before $T$, with $\tau = 0$ if there are no jumps.

This $\tau$ is not a stopping time, since knowing it needs a look at the future of the path, so we condition on all the jumps and on the path up to $\tau$. After $\tau$ there are no more jumps and the process is a Brownian motion with a modified drift, which means a Black and Scholes model with dividend yield $d + \lambda \kappa$. Hence the inner expectation is known in closed form, even if the process has jumps:
$$V_0 = \mathbb{E}(e^{-r \tau} C_{\lambda \kappa}(S_{\tau},T - \tau))$$
where $C_{\lambda \kappa}$ is the price of the previous section with $d$ replaced by $d + \lambda \kappa$. The discount is inside the expectation since $\tau$ is different for every path. This is what `ConditionalSaltando` computes.
```julia:conditional_saltando_derivatives
#hideall
@variables λ κ
with_jumps(delta) = substitute(delta, Dict(d => d + λ * κ))
println("\$\$\\partial_{S} C_{\\lambda \\kappa,\\text{call}} = ", latex_delta(with_jumps(delta_call)), "\$\$")
println("\$\$\\partial_{S} C_{\\lambda \\kappa,\\text{binary}} = ", latex_delta(with_jumps(delta_binary)), "\$\$")
```
With $h = T - \tau$, the same derivatives become:
\textoutput{conditional_saltando_derivatives}
On the same paths, this is the exact limit of Vibrato Saltando for a perfect inner integration. It has the same limits of the fixed $\tau$ case, it works only for European and binary options, and the same ones of Vibrato Saltando: $\lambda$ and $T$ can't be differentiated, since the last jump time depends on them.

For Variance Gamma and Normal Inverse Gaussian, `ConditionalSaltando` conditions on the subordinator $G_T$ instead (see Vibrato Saltando below): the inner price is Black and Scholes with volatility $\sigma \sqrt{G_T / T}$ and the rest of the drift as a dividend yield. There the price is a one dimensional integral over $G_T$, the formula of \cite{madan1998} for Variance Gamma, so the estimator is mostly a reference for `VibratoSaltando`.
\end{section}
\begin{section}{title="Vibrato Montecarlo Method"}
The pathwise method needs a smooth payoff, while the LRM has weights whose variance grows with the order of the sensitivity. Vibrato uses the best of each: it differentiates the path pathwise up to a time $\tau < T$, and applies the LRM only to what is left. The method was introduced in \cite{giles2009}, and its combination with automatic differentiation for high order sensitivities is studied in \cite{pages2018}.

Let's fix $\tau < T$, let's call $h = T - \tau$ and let's assume that the underlying process is a Levy process. Since:
$$V_0=\mathbb{E}(\Phi(X_T))=\mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))$$
$$\partial^{n}_{\theta} V_0=\mathbb{E}(\partial^{n}_{\theta} \mathbb{E}(\Phi(X_T)|X_{\tau}))$$
The inner expectation is smooth in $\theta$ even if $\Phi$ is not, since it is an integral of $\Phi$ against a smooth density. This is why the outer expectation can be differentiated. Applying the LRM to the inner expectation:
$$\partial^{n}_{\theta} V_0=\mathbb{E}(\mathbb{E}(\Phi(X_{T}) B_n(\{\partial^{i}_{\theta} \log{f_{X_T|X_{\tau}}(X_T,\theta)}\}_{i = 1 : n})|X_{\tau}))$$

As the traditional LRM, this needs the transition density. It is known for Ito processes, where the last step is Gaussian. For more general Levy processes it has to be inverted numerically from the characteristic function.

\begin{expand}{Gaussian last step}
If the last step is Gaussian, $X_T|X_{\tau} \sim \mathcal{N}(\mu, s^2)$, and $\mu$ and $s$ carry all the dependence on the parameters, with:
$$Z = \frac{X_T - \mu}{s}$$
then the scores with respect to $\mu$ and $s$, generated from `log_density`, are:
```julia:gaussian_scores
#hideall
@variables μ s
println("\$\$\\partial_{\\mu} \\log f = ", latex_weight(μ, s, μ, 1, x, Z), "\$\$")
println("\$\$\\partial_{s} \\log f = ", latex_weight(μ, s, s, 1, x, Z), "\$\$")
```
\textoutput{gaussian_scores}
and the score with respect to a generic parameter follows from the chain rule:
$$\partial_{\theta} \log{f} = \partial_{\mu} \log{f}\, \partial_{\theta} \mu + \partial_{s} \log{f}\, \partial_{\theta} s$$
In case of the Black and Scholes model:
$$\mu = X_{\tau} + (r-d-\frac{\sigma^2}{2}) h, \qquad s = \sigma \sqrt{h}$$
Note that $\mu$ depends on the parameters also through the path: this is the pathwise part of the method.
\end{expand}

\begin{expand}{Why automatic differentiation works}
The weights are never written by hand in the package. Let's call $\mathrm{v}(a)$ the value of a number $a$ without its derivatives. Then the sample points carry no derivative:
$$x(z) = \mathrm{v}(\mu) + \mathrm{v}(s)\, z$$
so the payoff is never differentiated, and all of the parameters enter through the density, in the factor:
$$m(g) = e^{g - \mathrm{v}(g)}$$
which has value one and, by Faa di Bruno's formula, derivatives $\partial^{n}_{\theta} m(g) = B_n(\{\partial^{i}_{\theta} g\}_{i = 1 : n})$. Hence evaluating at $g = \log{f(x(z);\mu,s)}$ gives the weights of every order out of any automatic differentiation number. The estimator on $N$ paths is, with $\mu_j$ and $s_j$ the values of the path $j$ and $c_j = \Phi(\mathrm{v}(\mu_j))$ the payoff at the centre of the density (see Payoff Centering below):
$$\hat{V}_0 = \frac{e^{-r T}}{N} \sum_{j=1}^{N} \int_{\mathbb{R}} \left((\Phi(x_j(z)) - c_j)\, m(\log{f(x_j(z);\mu_j,s_j)}) + c_j\right) \varphi(z) dz$$
where $\varphi$ is the standard normal density and $x_j(z) = \mathrm{v}(\mu_j) + \mathrm{v}(s_j) z$. The value is the plain Monte Carlo price, and its derivatives are the vibrato estimators.
\end{expand}

\begin{expand}{Weights for Black and Scholes}
```julia:vibrato_weights
#hideall
@variables x S_0 [latexwrapper = identity] r d σ T δ Z W_τ [latexwrapper = identity]
mu_bs = log(S_0) + (r - d - σ^2 / 2) * T + σ * W_τ
s_bs = σ * δ
println("\$\$w_{\\Delta} = ", latex_weight(mu_bs, s_bs, S_0, 1, x, Z), "\$\$")
println("\$\$w_{\\Gamma} = ", latex_weight(mu_bs, s_bs, S_0, 2, x, Z), "\$\$")
println("\$\$w_{\\mathcal{V}} = ", latex_weight(mu_bs, s_bs, σ, 1, x, Z), "\$\$")
```
Writing the path as $X_{\tau} = \log(S_0) + (r-d-\frac{\sigma^2}{2}) \tau + \sigma W_{\tau}$ and $\delta = \sqrt{h}$, we generate the weights with [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), starting from the same `log_density` used by the package and evaluating them at $X_T = \mu + \sigma \delta Z$:
\textoutput{vibrato_weights}
Each sensitivity is then, with $w$ the corresponding weight:
$$\partial^{n}_{\theta} V_0 = e^{-r T}\, \mathbb{E}(\mathbb{E}((\Phi(\mu + \sigma \delta Z) - \Phi(\mu))\, w|X_{\tau}))$$
For $\tau = 0$ these are the LRM weights of the previous section. The smaller $\delta$ is, the larger the weights are, and only the centering of the payoff keeps their variance under control. In $w_{\mathcal{V}}$ the term with $W_{\tau}$ comes from the path, while $Z^2 - 1$ comes from the last step alone.
\end{expand}

\end{section}

\begin{section}{title="Vibrato Saltando Montecarlo"}
Vibrato needs the law of the last step, and for a Levy process that is only available numerically. For Variance Gamma it can even be singular on a short step.
Vibrato Saltando \cite{scaramuzzino2017} chooses what to condition on so that the inner law is exactly Gaussian, for any number of time steps.

The method was introduced in the Master Thesis \cite{scaramuzzino2017}, Section 7.6, for finite activity Levy processes. "Vibrato" comes from the method of Giles \cite{giles2009}; "Saltando" (Italian for "jumping") refers both to the jumps (*salti*) of the processes it was made for and to a bowing technique of string players: to change note across the instrument, instead of keeping the bow on the strings they lift it and let it bounce on the next string.

\begin{expand}{Finite activity processes}
Let's assume that $X_t$ is a finite activity Levy process, like Kou or Merton:
$$X_t = \log(S_0) + (r-d-\frac{\sigma^2}{2}-\lambda \kappa) t + \sigma W_t + \sum_{i=1}^{N_t} Y_i, \qquad \kappa = \mathbb{E}(e^{Y_1} - 1)$$
and let's call $\tau$ the time of the last jump before $T$, with $\tau = 0$ if there are no jumps.
This $\tau$ is random and it is not a stopping time, since knowing it needs a look at the future of the path. So we condition on all the jumps and on $W_{\tau}$: the jumps are independent of $W$, and after $\tau$ there is nothing left but the diffusion. Calling $J = \sum_{i=1}^{N_T} Y_i$:
$$X_T|\tau, J, W_{\tau} \sim \mathcal{N}(\mu, s^2)$$
$$\mu = X_{\tau} + (r-d-\frac{\sigma^2}{2}-\lambda \kappa)(T - \tau), \qquad s = \sigma \sqrt{T - \tau}$$
The method is the one of the previous section, with a step that changes from path to path: the jumps go in the pathwise part, and the LRM is applied to the Gaussian diffusion after the last jump. Hence the payoff is integrated exactly, and $Z$ has no discretization error. On a path with no jumps $\tau = 0$, and Vibrato Saltando is the plain LRM over the whole maturity \cite{scaramuzzino2017}.
```julia:saltando_weights
#hideall
@variables J λ κ
mu_saltando = log(S_0) + (r - d - σ^2 / 2 - λ * κ) * T + σ * W_τ + J
println("\$\$w_{\\Delta} = ", latex_weight(mu_saltando, s_bs, S_0, 1, x, Z), "\$\$")
println("\$\$w_{\\Gamma} = ", latex_weight(mu_saltando, s_bs, S_0, 2, x, Z), "\$\$")
println("\$\$w_{\\mathcal{V}} = ", latex_weight(mu_saltando, s_bs, σ, 1, x, Z), "\$\$")
```
Generating the weights as before, with $\delta = \sqrt{T - \tau}$:
\textoutput{saltando_weights}
They are the same weights of Black and Scholes, since neither $J$ nor $\lambda \kappa$ depend on $S_0$ or $\sigma$. What changes is $\delta$, which is now different for every path.

Jumps are applied at the next point of the time grid, so the simulation returns the grid point before the step containing the last jump. This keeps the last step exact and at least one time step long: otherwise $T - \tau$ could be zero, with infinite weights.

The last jump time depends on $\lambda$ and the length of the last step on $T$, while the weights above consider them fixed: differentiating with respect to $\lambda$ or $T$ is not supported.
\end{expand}

\begin{expand}{Infinite activity processes}
Infinite activity processes like Variance Gamma \cite{madan1998} and Normal Inverse Gaussian \cite{barndorff1997} have no last jump, but they can be written as a Brownian motion with drift evaluated at a random time $G_t$ (the subordinator):
$$X_t = \log(S_0) + (r - d - \psi) t + \beta G_t + \sigma W_{G_t}$$
where $\beta$ and $\sigma$ are model parameters (called `θ` and `σ` in FinancialMonteCarlo), and $\psi$ is the drift correction such that $\mathbb{E}(S_T) = S_0 e^{(r-d)T}$.
Conditioning on the subordinator, the inner law is exactly Gaussian:
$$X_T|G_T \sim \mathcal{N}(\mu, s^2)$$
$$\mu = \log(S_0) + (r - d - \psi) T + \beta G_T, \qquad s = \sigma \sqrt{G_T}$$
Hence only $G_T$ has to be simulated (Gamma for Variance Gamma, Inverse Gaussian for Normal Inverse Gaussian), and the LRM applies to the Gaussian inner expectation with no time discretization.
```julia:subordinated_weights
#hideall
@variables β ψ
mu_sub = log(S_0) + (r - d - ψ) * T + β * δ^2
println("\$\$w_{\\Delta} = ", latex_weight(mu_sub, s_bs, S_0, 1, x, Z), "\$\$")
println("\$\$w_{\\Gamma} = ", latex_weight(mu_sub, s_bs, S_0, 2, x, Z), "\$\$")
```
The delta and gamma weights are once more the same, now with $\delta = \sqrt{G_T}$:
\textoutput{subordinated_weights}
For the vega, $\psi$ depends on $\sigma$ too: there the automatic differentiation of the factor $m$ takes care of it.

This is the extension to more general models that the thesis left as future work \cite{scaramuzzino2017}. Unlike $\lambda$ for finite activity processes, the parameters of the law of $G_T$ (its variance rate $\kappa$, and $T$) can be differentiated: $G_T$ is drawn from their values, and they enter through its density, in the same factor $m$ applied to the outer expectation:
$$V_0 = \mathbb{E}(\mathbb{E}(e^{-rT} \Phi(X_T)|G_T)\, m(\log{f_{G_T}(G_T;\kappa,T)}))$$
The draws of a Gamma or Inverse Gaussian sampler are not smooth in $\kappa$, while the density is. The noise of these derivatives is reduced by subtracting the sample mean of the inner expectations, as in Payoff Centering below.
\end{expand}

### Summary
All of the methods condition on something, and integrate the last part against its law. They differ in what they condition on, which gives the law of what is left ($\mu$ is its mean), and in how the inner expectation is computed:

| Process | Function | Conditioning on | Law of $X_T$ | Inner expectation |
|:--------|:---------|:----------------|:------------------------------------|:------------------|
| Black and Scholes | `Conditional` | path up to $T - h$ | $\mathcal{N}(\mu, \sigma^2 h)$ | closed form |
| Black and Scholes | `Vibrato` | path up to $T - h$ | $\mathcal{N}(\mu, \sigma^2 h)$ | integration |
| Heston \cite{heston1993} | `Vibrato` | path and variance $v$ up to $T - h$ | $\mathcal{N}(\mu, v_{T-h} h)$ | integration |
| Levy | `Vibrato` | path up to $T - h$ | FFT density | integration |
| Kou, Merton | `ConditionalSaltando` | jumps and path up to the last jump $\tau$ | $\mathcal{N}(\mu, \sigma^2 (T - \tau))$ | closed form |
| Kou, Merton | `VibratoSaltando` | jumps and path up to the last jump $\tau$ | $\mathcal{N}(\mu, \sigma^2 (T - \tau))$ | integration |
| Variance Gamma, Normal Inverse Gaussian | `ConditionalSaltando` | subordinator $G_T$ | $\mathcal{N}(\mu, \sigma^2 G_T)$ | closed form |
| Variance Gamma, Normal Inverse Gaussian | `VibratoSaltando` | subordinator $G_T$ | $\mathcal{N}(\mu, \sigma^2 G_T)$ | integration |


\end{section}

\begin{section}{title="Numerical Refinements"}
The previous sections are exact in theory; the following refinements make the inner expectation accurate in practice.
In all of them the inner law is Gaussian, $X_T | X_{\tau} \sim \mathcal{N}(\mu, s^2)$, and the inner expectation is computed on the standardized variable $Z = (X_T - \mu)/s$.

\begin{expand}{Payoff Centering}
For every $\theta$ the density integrates to one:
$$\int_{\mathbb{R}} f(x,\theta) dx = 1$$
hence all of its derivatives vanish and, for $n \geq 1$:
$$\mathbb{E}(B_n(\{\partial^{i}_{\theta} \log{f(X_T,\theta)}\}_{i = 1 : n})) = 0$$
Therefore any constant $c$ not depending on $\theta$ can be subtracted from the payoff without changing the sensitivities:
$$\partial^{n}_{\theta} V_0 = \mathbb{E}((\Phi(X_T) - c) B_n(\{\partial^{i}_{\theta} \log{f(X_T,\theta)}\}_{i = 1 : n}))$$
We take $c = \Phi(\mu)$, the payoff at the centre of the inner density, and add it back to the price.
The price is unchanged, while the score terms are now multiplied by $\Phi(\mu + s Z) - \Phi(\mu) \approx \Phi'(\mu) s Z$ instead of $\Phi(\mu + s Z)$.
In vibrato $s = \sigma \sqrt{T - \tau}$ is small and the score terms grow like $s^{-n}$, so without centering the variance of the delta estimator grows like $s^{-2}$; with centering it stays bounded.
Together with antithetic pairs $\pm Z$, the second order term becomes $\frac{1}{2}(\Phi(\mu + s Z) + \Phi(\mu - s Z)) - \Phi(\mu) \approx \frac{1}{2}\Phi''(\mu) s^2 Z^2$, which compensates the $s^{-2}$ of the gamma score.
In practice this reduces the noise of the gamma estimated with inner Monte Carlo by a factor 10 to 80.
\end{expand}

\begin{expand}{Independent Inner Samples}
The inner expectation is estimated with $M$ samples $Z_1, \dots, Z_M$ for each of the $N$ outer paths.
If the same samples were used for every outer path, the inner error would be the same for all paths and would not average out: the estimator would converge, as $N \to \infty$, to a biased value with error of order $M^{-1/2}$ (deltas 8-15% off with $M = 200$).
Drawing new inner samples for every outer path makes the inner errors independent, and the total error decreases with both $N$ and $M$.
With quasi Monte Carlo inner samples, the same Sobol points are randomly shifted modulo 1 for every outer path (Cranley-Patterson rotation \cite{cranley1976}), which keeps each inner estimate unbiased.
\end{expand}

\begin{expand}{Quadrature}
Alternatively, the inner expectation can be computed as a deterministic integral on a uniform grid in $Z$.
The grid must reach far enough in the tails: for a call the payoff grows exponentially in the upper tail, and truncating at $\pm 5$ standard deviations loses enough of the second order score terms to bias the volga by 3% to 45% on Variance Gamma. The default grid spans $\pm 7$ standard deviations.
\end{expand}

\end{section}

\begin{section}{title="American Options"}
An American option, exercisable at the dates $t_k = k\,h$, $k = 1, \dots, N$, $h = T/N$, has the price
$$V_0 = \sup_{\tau} \mathbb{E}(e^{-r \tau} \Phi(X_{\tau}))$$
The exercise policy is estimated with the Longstaff-Schwartz algorithm \cite{longstaff2001}: at each date the continuation value is regressed on $1, S, S^2$ over the in the money paths, and the option is exercised where the payoff exceeds it. The regression uses values only, so the policy carries no derivative.

\begin{expand}{Why freezing the exercise is not enough}
The usual approach \cite{pages2018} freezes the stopping time of each path and differentiates the cash flows as for a European contract. Since the policy is optimal, the first order effect of the moving exercise boundary vanishes and the first order sensitivities are correct. The second order ones are not: the boundary contributes a term proportional to the jump of the slope of the value at the boundary, and a frozen stopping time removes it. For a put, the frozen gamma only counts the paths that reach maturity near the strike, so it is even below the European gamma (-23% to -61% on the Longstaff-Schwartz put with $K = 40$).
\end{expand}

\begin{expand}{Vibrato over every exercise date}
A policy frozen in price space, i.e. an exercise region $E_k$ of prices at each date, does not have this problem for the sensitivities with respect to $S_0$: the optimal exercise boundary does not depend on $S_0$, so the policy is optimal for every $S_0$, and the error of the regressed boundary enters every derivative with respect to $S_0$ only at second order.
The discontinuity of the exercise decision is then handled as the payoff at maturity: by vibrato. The paths are simulated at $s_0 = 0$ and $s_k = (k - \frac{1}{2}) h$, so that every exercise date is the midpoint of $[s_k, s_{k+1}]$ and, given the ends, the log price at $t_k$ is a Brownian bridge:
$$X_{t_k} | X_{s_k}, X_{s_{k+1}} \sim \mathcal{N}(\frac{X_{s_k} + X_{s_{k+1}}}{2}, \frac{\sigma^2 h}{4})$$
The bridge values are conditionally independent, and the future of a path depends on $X_{t_k}$ only through the decision at $t_k$. The value of each path is computed backwards:
$$W_N = \mathbb{E}(e^{-r T} \Phi(X_T) | X_{s_N})$$
$$W_k = \mathbb{E}(1_{X_{t_k} \in E_k} e^{-r t_k} \Phi(X_{t_k}) + 1_{X_{t_k} \notin E_k} W_{k+1} | X_{s_k}, X_{s_{k+1}})$$
and $V_0 = \mathbb{E}(W_1)$. Every conditional expectation is computed with the LRM over the bridge density, so the exercise indicator and the payoff are never differentiated: the estimator is smooth in all the parameters and unbiased for the chosen policy. The price of this exactness is variance: at the exercise boundary the integrand jumps by the gap between the exercise value and the realized continuation, and the weights of a bridge of standard deviation $\sigma \sqrt{h} / 2$ amplify it, more and more as the number of exercise dates grows.

The bridge is Gaussian only for Black-Scholes dynamics. For Lévy and Heston models the exercise decision for $t_k$ is taken at $t_{k-1}$, the payoff is replaced by its vibrato expectation over the step, and the decisions are frozen: their first order sensitivities are accurate, the second order ones are biased.
\end{expand}

\end{section}

\begin{section}{title="References"}
* \biblabel{barndorff1997}{Barndorff-Nielsen (1997)} O. E. Barndorff-Nielsen, *Processes of normal inverse Gaussian type*, Finance and Stochastics 2(1), 41-68, 1997.
* \biblabel{blackscholes1973}{Black, Scholes (1973)} F. Black, M. Scholes, *The pricing of options and corporate liabilities*, Journal of Political Economy 81(3), 637-654, 1973.
* \biblabel{broadie1996}{Broadie, Glasserman (1996)} M. Broadie, P. Glasserman, *Estimating security price derivatives using simulation*, Management Science 42(2), 269-285, 1996. [doi](https://doi.org/10.1287/mnsc.42.2.269)
* \biblabel{cranley1976}{Cranley, Patterson (1976)} R. Cranley, T. N. L. Patterson, *Randomization of number theoretic methods for multiple integration*, SIAM Journal on Numerical Analysis 13(6), 904-914, 1976. [doi](https://doi.org/10.1137/0713071)
* \biblabel{giles2009}{Giles (2009)} M. B. Giles, *Vibrato Monte Carlo sensitivities*, in P. L'Ecuyer, A. B. Owen (eds.), Monte Carlo and Quasi-Monte Carlo Methods 2008, 369-382, Springer, 2009. [doi](https://doi.org/10.1007/978-3-642-04107-5_23)
* \biblabel{glasserman2004}{Glasserman (2004)} P. Glasserman, *Monte Carlo Methods in Financial Engineering*, Springer, 2004.
* \biblabel{heston1993}{Heston (1993)} S. L. Heston, *A closed-form solution for options with stochastic volatility with applications to bond and currency options*, The Review of Financial Studies 6(2), 327-343, 1993.
* \biblabel{kou2002}{Kou (2002)} S. G. Kou, *A jump-diffusion model for option pricing*, Management Science 48(8), 1086-1101, 2002.
* \biblabel{longstaff2001}{Longstaff, Schwartz (2001)} F. A. Longstaff, E. S. Schwartz, *Valuing American options by simulation: a simple least-squares approach*, The Review of Financial Studies 14(1), 113-147, 2001.
* \biblabel{madan1998}{Madan, Carr, Chang (1998)} D. B. Madan, P. P. Carr, E. C. Chang, *The variance gamma process and option pricing*, European Finance Review 2(1), 79-105, 1998.
* \biblabel{merton1976}{Merton (1976)} R. C. Merton, *Option pricing when underlying stock returns are discontinuous*, Journal of Financial Economics 3(1-2), 125-144, 1976.
* \biblabel{pages2018}{Pagès, Pironneau, Sall (2018)} G. Pagès, O. Pironneau, G. Sall, *Vibrato and automatic differentiation for high order derivatives and sensitivities of financial options*, Journal of Computational Finance 22(2), 2018. [arXiv](https://arxiv.org/abs/1606.06143)
* \biblabel{scaramuzzino2017}{Scaramuzzino (2017)} N. Scaramuzzino, *Differenziazione automatica: teoria ed applicazioni in finanza*, Master Thesis, Politecnico di Milano, supervisor C. Sgarra, academic year 2016/2017. Vibrato Saltando is in Section 7.6. [handle](https://hdl.handle.net/10589/135858)
\end{section}
