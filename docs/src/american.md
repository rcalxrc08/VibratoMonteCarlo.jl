# American options

An American option, exercisable at the dates ``t_k = k\,h``, ``k = 1, \dots, N``, ``h = T/N``, has the price

```math
V_0 = \sup_{\tau} \mathbb{E}(e^{-r \tau} \Phi(X_{\tau}))
```

The exercise policy is estimated with the Longstaff-Schwartz algorithm [longstaff2001](@cite): at each date the continuation value is regressed on ``1, S, S^2`` over the in the money paths, and the option is exercised where the payoff exceeds it. The regression uses values only, so the policy carries no derivative.

## Why freezing the exercise is not enough
The usual approach [pages2018](@cite) freezes the stopping time of each path and differentiates the cash flows as for a European contract. Since the policy is optimal, the first order effect of the moving exercise boundary vanishes and the first order sensitivities are correct. The second order ones are not: the boundary contributes a term proportional to the jump of the slope of the value at the boundary, and a frozen stopping time removes it. For a put, the frozen gamma only counts the paths that reach maturity near the strike, so it is even below the European gamma (-23% to -61% on the Longstaff-Schwartz put with ``K = 40``).

## Vibrato over every exercise date
A policy frozen in price space, i.e. an exercise region ``E_k`` of prices at each date, does not have this problem for the sensitivities with respect to ``S_0``: the optimal exercise boundary does not depend on ``S_0``, so the policy is optimal for every ``S_0``, and the error of the regressed boundary enters every derivative with respect to ``S_0`` only at second order.
The discontinuity of the exercise decision is then handled as the payoff at maturity: by vibrato. The paths are simulated at ``s_0 = 0`` and ``s_k = (k - \frac{1}{2}) h``, so that every exercise date is the midpoint of ``[s_k, s_{k+1}]`` and, given the ends, the log price at ``t_k`` is a Brownian bridge:

```math
X_{t_k} | X_{s_k}, X_{s_{k+1}} \sim \mathcal{N}(\frac{X_{s_k} + X_{s_{k+1}}}{2}, \frac{\sigma^2 h}{4})
```

The bridge values are conditionally independent, and the future of a path depends on ``X_{t_k}`` only through the decision at ``t_k``. The value of each path is computed backwards:

```math
W_N = \mathbb{E}(e^{-r T} \Phi(X_T) | X_{s_N})
```

```math
W_k = \mathbb{E}(1_{X_{t_k} \in E_k} e^{-r t_k} \Phi(X_{t_k}) + 1_{X_{t_k} \notin E_k} W_{k+1} | X_{s_k}, X_{s_{k+1}})
```

and ``V_0 = \mathbb{E}(W_1)``. Every conditional expectation is computed with the LRM over the bridge density, so the exercise indicator and the payoff are never differentiated: the estimator is smooth in all the parameters and unbiased for the chosen policy. The price of this exactness is variance: at the exercise boundary the integrand jumps by the gap between the exercise value and the realized continuation, and the weights of a bridge of standard deviation ``\sigma \sqrt{h} / 2`` amplify it, more and more as the number of exercise dates grows.

The bridge is Gaussian only for Black-Scholes dynamics. For Lévy and Heston models the exercise decision for ``t_k`` is taken at ``t_{k-1}``, the payoff is replaced by its vibrato expectation over the step, and the decisions are frozen: their first order sensitivities are accurate, the second order ones are biased.
