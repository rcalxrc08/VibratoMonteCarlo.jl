<!-- =============================
     ABOUT
    ============================== -->
\begin{section}{title="About this Package", name="About"}
\lead{VibratoMonteCarlo.jl is a Julia package containing some useful Financial functions for sensitivity computation.}
VibratoMonteCarlo.jl is built on top of [FinancialMonteCarlo.jl](https://github.com/JuliaDiff/DualNumbers.jl) and [FinancialFFT.jl](https://github.com/JuliaDiff/DualNumbers.jl).
Standard montecarlo methods lacks of differentiability, which makes automatic differentiation useless.
The main aim of this package is to provide a feasible way to compute sensitivities of any order for various types of payoffs using montecarlo methods.

The implementation for Vibrato Montecarlo is heavily based on the Master Thesis of Nicola Scaramuzzino (that is me) and the subsequent development.

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

If we develop the integral:

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
$$B_0 =   1 $$
$$B_1(x_1) =   x_1 $$
$$B_2(x_1,x_2) =   x_1^2 + x_2 $$
$$B_3(x_1,x_2,x_3) =   x_1^3 + 3x_1 x_2 + x_3 $$
$$B_4(x_1,x_2,x_3,x_4) =   x_1^4 + 6 x_1^2 x_2 + 4 x_1 x_3 + 3 x_2^2 + x_4 $$
$$B_5(x_1,x_2,x_3,x_4,x_5) =   x_1^5 + 10 x_2 x_1^3 + 15 x_2^2 x_1 + 10 x_3 x_1^2 + 10 x_3 x_2 + 5 x_4 x_1 + x_5 $$
$$B_6(x_1,x_2,x_3,x_4,x_5,x_6) =   x_1^6 + 15 x_2 x_1^4 + 20 x_3 x_1^3 + 45 x_2^2 x_1^2 + 15 x_2^3 + 60 x_3 x_2 x_1 + 15 x_4 x_1^2 + 10 x_3^2 + 15 x_4 x_2 + 6 x_5 x_1 + x_6 $$
Hence:

$$\partial^{n}_{\theta} f(x,\theta) = f(x,\theta) B_n(\{\partial^{i}_{\theta} \log{f(x,\theta)}\}_{i = 1 : n})  $$

By applying this to the original problem:

$$\partial^{n}_{\theta} V_0=\partial^{n}_{\theta} \int_{\mathbb{R}}  \Phi(x) f(x,\theta) dx = \int_{\mathbb{R}}  \Phi(x) \partial^{n}_{\theta} f(x,\theta) dx$$

By applying the previous result:

$$\partial^{n}_{\theta} V_0= \int_{\mathbb{R}}  \Phi(x) f(x,\theta) B_n(\{\partial^{i}_{\theta} \log{f(x,\theta)}\}_{i = 1 : n}) dx$$

or, in terms of expectation:


$$\partial^{n}_{\theta} V_0 = \mathbb{E}(\Phi(X_T) B_n(\{\partial^{i}_{\theta} \log{f(X_T,\theta)}\}_{i = 1 : n}))$$




### Likelihood Ratio Method for Black and Scholes
In case of trivial models the density is known in analytic form, otherwise one need to compute it numerically from the characteristic function.
In case of the Black and Scholes Model, and a European Option we have the following:

$$X_t=\log(S_0)+(r-d-\dfrac{\sigma^2}{2}) t + \sigma W_t$$

$$\log(f(x,S_0,r,d,\sigma,T)) = -\dfrac{1}{2}((\dfrac{x-(\log(S_0)+(r-d-\dfrac{\sigma^2}{2}) T)}{\sigma \sqrt{T}})^2 + \log(2 \pi)) - \log(\sigma\sqrt{T})$$

From this we can compute the various partial derivatives:

$$ \partial_{S_0} \log(f(x,S_0,r,d,\sigma,T)) = \frac{x - \log\left( S_0 \right) - \left(r - d - \frac{\sigma^{2}}{2} \right) T}{\sigma^{2} T S_0}$$
$$ \partial_{r} \log(f(x,S_0,r,d,\sigma,T)) = \frac{ x - \log\left( S_0 \right) - T \left(  r - d - \frac{\sigma^{2}}{2} \right) }{\sigma^{2}}$$
$$ \partial_{d} \log(f(x,S_0,r,d,\sigma,T)) = - \partial_{r} \log(f(x,S_0,r,d,\sigma,T))$$
$$ \partial_{\sigma} \log(f(x,S_0,r,d,\sigma,T)) =\frac{ - \left( x - \log\left( S_0 \right) - T \left(  r - d - \frac{\sigma^{2}}{2} \right) \right) \left( \frac{T}{\sqrt{T}} - \frac{x - \log\left( S_0 \right) - T \left(  r - d - \frac{\sigma^{2}}{2} \right)}{\sigma^{2} T} \sqrt{T} \right)}{\sqrt{T} \sigma} + \frac{-1}{\sigma}$$
$$ \partial_{T} \log(f(x,S_0,r,d,\sigma,T)) =\frac{ - \left( x - \log\left( S_0 \right) - T \left(  r - d - \frac{\sigma^{2}}{2} \right) \right) \left( \frac{d - r + \frac{\sigma^{2}}{2}}{\sqrt{T} \sigma} + \frac{ - \frac{x - \log\left( S_0 \right) - T \left(  r - d - \frac{\sigma^{2}}{2} \right)}{\sigma^{2} T} \sigma}{2 \sqrt{T}} \right)}{\sqrt{T} \sigma} + \frac{-1}{2 T}$$

Hence, a delta for a generic european option can be computed as the following expectation:

$$\Delta = \partial_{S_0} V_0 = \mathbb{E}(\Phi(X_T) \frac{X_T - \log\left( S_0 \right) - \left(r - d - \frac{\sigma^{2}}{2} \right) T}{\sigma^{2} T S_0} )$$

or, in terms of $W_T$:

$$\Delta = \partial_{S_0} V_0 = \mathbb{E}(\Phi(\log(S_0)+(r-d-\dfrac{\sigma^2}{2}) T + \sigma W_T) \frac{W_T}{\sigma T S_0} )$$

Let's assume now that $\Phi$ is a forward i.e.:
$$\Phi(x)=e^x$$

Since $\Delta = e^{(r-d)T}$:
$$\mathbb{E}(e^{-\frac{\sigma^2}{2} T + \sigma W_T}\,W_T) = \sigma T$$

## Spectral Theorem for Black and Scholes

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

### Application for sixth order derivative
```julia:pure_derivative
#hideall
using Symbolics,Latexify,VibratoMonteCarlo
function generate_and_adjust(der_sigma,mu,sigma,S_0,r,d,σ,T)
    first_dict=Dict(mu=>log(S_0)+(r-d-σ^2/2)*T,sigma=>σ*sqrt(T));
    new_der_sigma=substitute(der_sigma,first_dict)
	return new_der_sigma
end

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


\end{section}
\begin{section}{title="Conditional Expectation"}
Let's fix a positive number $\tau < T$ and let's assume now that the underlying process is a Levy process.\\
Since:
$$V_0=\mathbb{E}(\Phi(X_T))=\mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))$$
$$\partial_{\theta} V_0=\mathbb{E}(\Phi(X_T))=\partial_{\theta} \mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))=\mathbb{E}(\partial_{\theta} \mathbb{E}(\Phi(X_{T})|X_{\tau}))$$
\end{section}
\begin{section}{title="Conditional Saltando Expectation"}
Let's assume $\tau$ is a stopping time st $\tau < T$ a.s., and let's assume now that the underlying process is a Levy process.\\
Since:
$$V_0=\mathbb{E}(\Phi(X_T))=\mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))$$
$$\partial_{\theta} V_0=\mathbb{E}(\Phi(X_T))=\partial_{\theta} \mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))=\mathbb{E}(\partial_{\theta} \mathbb{E}(\Phi(X_{T})|X_{\tau}))$$

Furthemore, if we assume that $X_t$ is a finite activity Levy process, a smart choice of the law of the stopping time can benefit the computation of the inner expectation.
Indeed, if we take $\tau$ as the stopping time corresponding to the last jump of the process, we have that:
$$ X_t | X_{\tau}  \text{ is a Ito process }$$
Why is it useful?\\
Let's assume $X_t$ is a Kou process, then $X_t | X_{\tau}$ is a Brownian motion with a modified drift, hence we can have an analytic formula for the inner expectation.
\end{section}
\begin{section}{title="Vibrato Montecarlo Method"}
Let's fix a positive number $\tau < T$ and let's assume now that the underlying process is a Levy process.\\
Since:
$$V_0=\mathbb{E}(\Phi(X_T))=\mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))=\mathbb{E}(\mathbb{E}(\Phi(X_{\tau}+(X_{T}-X_{\tau}))|X_{\tau}))$$
$$\partial_{\theta} V_0=\mathbb{E}(\Phi(X_T))=\partial_{\theta} \mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))=\partial_{\theta} \mathbb{E}(\mathbb{E}(\Phi(X_{\tau}+(X_{T}-X_{\tau}))|X_{\tau}))$$
$$\partial_{\theta} V_0=\mathbb{E}(\partial_{\theta} \mathbb{E}(\Phi(X_{\tau}+(X_{T}-X_{\tau}))|X_{\tau}))$$

Now we apply the LRM method to the inner expectation, let's call $Y_{T-\tau}= X_{T}-X_{\tau}$, then we have:

$$\partial_{\theta} \mathbb{E}(\Phi(X_{T})|X_{\tau}) = \mathbb{E}(\Phi(X_{T})\partial_{\theta}  \log(f_{X_T-X_{\tau}}(X_T - X_{\tau},\theta))|X_{\tau}) $$

Hence
$$\partial_{\theta} V_0=\mathbb{E}(\mathbb{E}(\Phi(X_{T})\partial_{\theta}  \log(f_{X_T-X_{\tau}}(X_T - X_{\tau},\theta))|X_{\tau}))$$

or for a generic order n:

$$\partial^{n}_{\theta} V_0=\mathbb{E}(\mathbb{E}(\Phi(X_{T}) B_n(\{\partial^{i}_{\theta} \log{f_{X_T-X_{\tau}}(X_T - X_{\tau},\theta)}\}_{i = 1 : n})))$$


As the traditional LR method, this is hard to apply in case the transition density is unknown. This can be easily solved in case of Ito processes where we can approximate the increment as a sum of powers of normal random variables.
In case of more general process like Levy, one needs to invert numerically the density from the characteristic function.

\end{section}

\begin{section}{title="Vibrato Saltando Montecarlo"}
Since what we stated above applies also to a generic stopping time $\tau$, we investigate some smart choice of the law of $\tau$ in order to exploit some property of the underlying process.
Let's assume now that $X_t$ is a Finite Activity Levy process and that $\tau$ is the stopping time corresponding to the last jump of the underlying.
We know already that $X_t|X_{\tau}$ is a Ito process.
On such process we are able to apply the standard LRM or VBM since we know the corresponding law of the increments (powers of gaussians).

\end{section}

\begin{section}{title="Numerical Refinements"}
The previous sections are exact in theory; the following refinements make the inner expectation accurate in practice.
In all of them the inner law is Gaussian, $X_T | X_{\tau} \sim \mathcal{N}(\mu, s^2)$, and the inner expectation is computed on the standardized variable $Z = (X_T - \mu)/s$.

### Payoff Centering
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

### Independent Inner Samples
The inner expectation is estimated with $M$ samples $Z_1, \dots, Z_M$ for each of the $N$ outer paths.
If the same samples were used for every outer path, the inner error would be the same for all paths and would not average out: the estimator would converge, as $N \to \infty$, to a biased value with error of order $M^{-1/2}$ (deltas 8-15% off with $M = 200$).
Drawing new inner samples for every outer path makes the inner errors independent, and the total error decreases with both $N$ and $M$.
With quasi Monte Carlo inner samples, the same Sobol points are randomly shifted modulo 1 for every outer path (Cranley-Patterson rotation), which keeps each inner estimate unbiased.

### Quadrature
Alternatively, the inner expectation can be computed as a deterministic integral on a uniform grid in $Z$.
The grid must reach far enough in the tails: for a call the payoff grows exponentially in the upper tail, and truncating at $\pm 5$ standard deviations loses enough of the second order score terms to bias the volga by 3% to 45% on Variance Gamma. The default grid spans $\pm 7$ standard deviations.

### Vibrato Saltando for Subordinated Processes
Infinite activity processes like Variance Gamma and Normal Inverse Gaussian have no last jump, but they can be written as a Brownian motion with drift evaluated at a random time $G_t$ (the subordinator):
$$X_t = \log(S_0) + (r - d - \psi) t + \beta G_t + \sigma W_{G_t}$$
where $\beta$ and $\sigma$ are model parameters (called `θ` and `σ` in FinancialMonteCarlo), and $\psi$ is the drift correction such that $\mathbb{E}(S_T) = S_0 e^{(r-d)T}$.
Conditioning on the subordinator instead of on a stopping time, the inner law is exactly Gaussian:
$$X_T | G_T \sim \mathcal{N}(\log(S_0) + (r - d - \psi) T + \beta G_T, \sigma^2 G_T)$$
Hence only $G_T$ has to be simulated (Gamma for Variance Gamma, Inverse Gaussian for Normal Inverse Gaussian), and the LRM applies to the Gaussian inner expectation with no time discretization.
\end{section}

\begin{section}{title="American Options"}
An American option, exercisable at the dates $t_k = k\,h$, $k = 1, \dots, N$, $h = T/N$, has the price
$$V_0 = \sup_{\tau} \mathbb{E}(e^{-r \tau} \Phi(X_{\tau}))$$
The exercise policy is estimated with the Longstaff-Schwartz algorithm: at each date the continuation value is regressed on $1, S, S^2$ over the in the money paths, and the option is exercised where the payoff exceeds it. The regression uses values only, so the policy carries no derivative.

### Why freezing the exercise is not enough
The usual approach freezes the stopping time of each path and differentiates the cash flows as for a European contract. Since the policy is optimal, the first order effect of the moving exercise boundary vanishes and the first order sensitivities are correct. The second order ones are not: the boundary contributes a term proportional to the jump of the slope of the value at the boundary, and a frozen stopping time removes it. For a put, the frozen gamma only counts the paths that reach maturity near the strike, so it is even below the European gamma (-23% to -61% on the Longstaff-Schwartz put with $K = 40$).

### Vibrato over every exercise date
A policy frozen in price space, i.e. an exercise region $E_k$ of prices at each date, does not have this problem for the sensitivities with respect to $S_0$: the optimal exercise boundary does not depend on $S_0$, so the policy is optimal for every $S_0$, and the error of the regressed boundary enters every derivative with respect to $S_0$ only at second order.
The discontinuity of the exercise decision is then handled as the payoff at maturity: by vibrato. The paths are simulated at $s_0 = 0$ and $s_k = (k - \frac{1}{2}) h$, so that every exercise date is the midpoint of $[s_k, s_{k+1}]$ and, given the ends, the log price at $t_k$ is a Brownian bridge:
$$X_{t_k} | X_{s_k}, X_{s_{k+1}} \sim \mathcal{N}(\frac{X_{s_k} + X_{s_{k+1}}}{2}, \frac{\sigma^2 h}{4})$$
The bridge values are conditionally independent, and the future of a path depends on $X_{t_k}$ only through the decision at $t_k$. The value of each path is computed backwards:
$$W_N = \mathbb{E}(e^{-r T} \Phi(X_T) | X_{s_N})$$
$$W_k = \mathbb{E}(1_{X_{t_k} \in E_k} e^{-r t_k} \Phi(X_{t_k}) + 1_{X_{t_k} \notin E_k} W_{k+1} | X_{s_k}, X_{s_{k+1}})$$
and $V_0 = \mathbb{E}(W_1)$. Every conditional expectation is computed with the LRM over the bridge density, so the exercise indicator and the payoff are never differentiated: the estimator is smooth in all the parameters and unbiased for the chosen policy. The price of this exactness is variance: at the exercise boundary the integrand jumps by the gap between the exercise value and the realized continuation, and the weights of a bridge of standard deviation $\sigma \sqrt{h} / 2$ amplify it, more and more as the number of exercise dates grows.

The bridge is Gaussian only for Black-Scholes dynamics. For Lévy and Heston models the exercise decision for $t_k$ is taken at $t_{k-1}$, the payoff is replaced by its vibrato expectation over the step, and the decisions are frozen: their first order sensitivities are accurate, the second order ones are biased.
\end{section}