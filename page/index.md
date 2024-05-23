<!-- =============================
     ABOUT
    ============================== -->
\begin{section}{title="About this Package", name="About"}
\lead{VibratoMonteCarlo.jl is a Julia package containing some useful Financial functions for sensitivity computation.}
VibratoMonteCarlo.jl is built on top of [FinancialMonteCarlo.jl](https://github.com/JuliaDiff/DualNumbers.jl) and [FinancialFFT.jl](https://github.com/JuliaDiff/DualNumbers.jl).
Standard montecarlo methods lacks of differentiability, which makes automatic differentiation useless.
The main aim of this package is to provide a feasible way to compute sensitivities of any order for various types of payoffs using montecarlo methods.
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

$$\partial_\theta V_0=\partial_\theta \mathbb{E}(\Phi(X_T)) = \partial_\theta \mathbb{E}(\dfrac{\partial \Phi(X_T)}{\partial X_T}\dfrac{\partial X_T}{\partial \theta})$$

Where the following is called tangent process:

$$Y_t=\dfrac{\partial X_t}{\partial \theta}$$

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

$$\partial_{\theta^n} e^{g(x,\theta)} = e^{g(x,\theta)} B_n(\{\partial_{\theta^i} g(x,\theta)\}_{i = 1 : n})  $$

Where $B_n$ is the complete exponential [Bell polynomial](https://en.wikipedia.org/wiki/Bell_polynomials) of n-th order.
$$B_0 =   1 $$
$$B_1(x_1) =   x_1 $$
$$B_2(x_1,x_2) =   x_1^2 + x_2 $$
$$B_3(x_1,x_2,x_3) =   x_1^3 + 3x_1 x_2 + x_3 $$
$$B_4(x_1,x_2,x_3,x_4) =   x_1^4 + 6 x_1^2 x_2 + 4 x_1 x_3 + 3 x_2^2 + x_4 $$
$$B_5(x_1,x_2,x_3,x_4,x_5) =   x_1^5 + 10 x_2 x_1^3 + 15 x_2^2 x_1 + 10 x_3 x_1^2 + 10 x_3 x_2 + 5 x_4 x_1 + x_5 $$
$$B_6(x_1,x_2,x_3,x_4,x_5,x_6) =   x_1^6 + 15 x_2 x_1^4 + 20 x_3 x_1^3 + 45 x_2^2 x_1^2 + 15 x_2^3 + 60 x_3 x_2 x_1 + 15 x_4 x_1^2 + 10 x_3^2 + 15 x_4 x_2 + 6 x_5 x_1 + x_6 $$
Hence:

$$\partial_{\theta^n} f(x,\theta) = f(x,\theta) B_n(\{\partial_{\theta^i} \log{f(x,\theta)}\}_{i = 1 : n})  $$

By applying this to the original problem:

$$\partial_{\theta^n} V_0=\partial_{\theta^n} \int_{\mathbb{R}}  \Phi(x) f(x,\theta) dx = \int_{\mathbb{R}}  \Phi(x) \partial_{\theta^n} f(x,\theta) dx$$

By applying the previous result:

$$\partial_{\theta^n} V_0= \int_{\mathbb{R}}  \Phi(x) f(x,\theta) B_n(\{\partial_{\theta^i} \log{f(x,\theta)}\}_{i = 1 : n}) dx$$

or, in terms of expectation:


$$\partial_{\theta^n} V_0 = \mathbb{E}(\Phi(X_T) B_n(\{\partial_{\theta^i} \log{f(X_T,\theta)}\}_{i = 1 : n}))$$



\end{section}


\begin{section}{title="Likelihood Ratio Method for Black and Scholes"}
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

\end{section}
\begin{section}{title="Spectral Theorem for Black and Scholes"}

Let's observe something, $\log(f(x,\theta))$ is a second order polymial in $x$:

$$\log{f(X_T,\theta)}= \sum_{i=0}^2 b_i(\theta) X_T^i = b_0(\theta) + b_1(\theta)\,X_T +  b_2(\theta)\,X_T^2$$

Then all al of the partial derivatives will be second order polynomials as well:

$$\partial_{\theta^i} \log{f(X_T,\theta)}= \partial_{\theta^i}b_0(\theta) + \partial_{\theta^i}b_1(\theta)\,X_T +  \partial_{\theta^i}b_2(\theta)\,X_T^2$$

And the Bell exponential polymial will be a polynomial of degree $2 n$ in terms of $X_T$:

$$B_n(\{\partial_{\theta^i} \log{f(X_T,\theta)}\}_{i = 1 : n}) = \sum_{i=0}^{2 n} a_i(\theta,n) X_T^i $$

And finally:

$$\partial_{\theta^n} V_0  = a_0(\theta,n) V_0 + \sum_{i=1}^{2 n} a_i(\theta,n) \mathbb{E}(\Phi(X_T) X_T^i) $$

\alert{**Spectral Theorem**: Given a European Option with payoff $\Phi$, an underlying varying according to Black and Scholes model, then we can express a generic $n$ order sensitivity as follows:
$$\partial_{\theta^n} V_0  = a_0(\theta,n) V_0 + \sum_{i=1}^{2 n} a_i(\theta,n) \mathbb{E}(\Phi(X_T) X_T^i) $$

Hence a generic $n$ order sensitivity can be express as a linear combination of the prices of $2 n + 1$ contracts.
}
To be noticed that the various coefficients $a_i(\theta,n)$ do not depend on the payoff, hence with $2n + 1$ prices we are able to compute all of the partial derivatives with respect to any parameter up to order $n$.
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

$$\partial_{\theta^n} V_0=\mathbb{E}(\mathbb{E}(\Phi(X_{T}) B_n(\{\partial_{\theta^i} \log{f_{X_T-X_{\tau}}(X_T - X_{\tau},\theta)}\}_{i = 1 : n})))$$


As the traditional LR method, this is hard to apply in case the transition density is unknown. This can be easily solved in case of Ito processes where we can approximate the increment as a sum of powers of normal random variables.
In case of more general process like Levy, one needs to invert numerically the density from the characteristic function.

\end{section}

\begin{section}{title="Vibrato Saltando Montecarlo"}
Since what we stated above applies also to a generic stopping time $\tau$, we investigate some smart choice of the law of $\tau$ in order to exploit some property of the underlying process.
Let's assume now that $X_t$ is a Finite Activity Levy process and that $\tau$ is the stopping time corresponding to the last jump of the underlying.
We know already that $X_t|X_{\tau}$ is a Ito process.
On such process we are able to apply the standard LRM or VBM since we know the corresponding law of the increments (powers of gaussians).

\end{section}