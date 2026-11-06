# Formulas of the documentation, generated with Symbolics at build time from VibratoMonteCarlo.log_density: no derived
# formula is typed by hand. Each function returns Markdown with display math, for an @eval block of the pages.
# δ = √h is a plain symbol, since Symbolics does not cancel sqrt(T).
using Symbolics, Latexify, VibratoMonteCarlo, Markdown

display_math(lines) = Markdown.parse(join(("```math\n$l\n```" for l in lines), "\n\n"))

# From the recurrence B_{n+1} = sum_{k=0}^n binomial(n, k) B_{n-k} x_{k+1}.
function generate_bell_polynomials(max_order)
    x_vars = Symbolics.variables(:x, 1:max_order)
    bell = Num[1]
    for n = 0:(max_order-1)
        push!(bell, expand(sum(binomial(n, k) * bell[n-k+1] * x_vars[k+1] for k = 0:n)))
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
function bell_polynomials()
    bell, x_vars = generate_bell_polynomials(6)
    return display_math(["B_0 = 1"; ["B_$n(" * join(("x_$i" for i = 1:n), ",") * ") = " * latex_bell(bell[n+1], x_vars[1:n]) for n = 1:6]])
end

# The score is evaluated at x = X_T before simplifying: what is left is a small polynomial in W_T,
# while simplifying the derivative in x first takes hours with SymbolicUtils 4.
function generate_faa_di_bruno(z_better, d_tgt_var, x, X_T)
    init_val = exp(z_better)
    score = expand_derivatives(d_tgt_var(init_val)) / init_val
    return simplify_fractions(expand(substitute(score, Dict(x => X_T))))
end
generate_faa_di_bruno_exp(z_better, tgt_var, order, x, X_T) = generate_faa_di_bruno(z_better, Differential(tgt_var)^order, x, X_T)
generate_faa_di_bruno_exp_mix(z_better, tgt_var, order, tgt_var2, order2, x, X_T) = generate_faa_di_bruno(z_better, Differential(tgt_var)^order * Differential(tgt_var2)^order2, x, X_T)

function generate_weight(mu, s, tgt_var, order, x, Z)
    z = VibratoMonteCarlo.log_density(x, mu, s)
    return generate_faa_di_bruno_exp(z, tgt_var, order, x, mu + s * Z)
end
latex_weight(weight) = latexify(weight, env = :raw, mult_symbol = "")
latex_weight(mu, s, tgt_var, order, x, Z) = latex_weight(generate_weight(mu, s, tgt_var, order, x, Z))

generate_and_adjust(der_sigma, mu, sigma, S_0, r, d, σ, T) = substitute(der_sigma, Dict(mu => log(S_0) + (r - d - σ^2 / 2) * T, sigma => σ * sqrt(T)))

function generate_latex_expression_base(S_0, r, d, σ, T, tgt_var, max_order)
    @variables mu, sigma, x
    z = VibratoMonteCarlo.log_density(x, mu, sigma)
    @variables W_T [latexwrapper = identity]
    X_T = log(S_0) + (r - d - σ^2 / 2) * T + σ * W_T
    z_better = generate_and_adjust(z, mu, sigma, S_0, r, d, σ, T)
    adjusted_final_der = generate_faa_di_bruno_exp(z_better, tgt_var, max_order, x, X_T)
    @variables V ∂ ∂θ [latexwrapper = identity]
    der_num = ∂^max_order * V / (∂θ^max_order)
    @variables Φ(X_T)
    ress_ = Φ * adjusted_final_der
    @variables 𝔼(ress_)
    return latexify(der_num ~ 𝔼, env = :raw, mult_symbol = "").s
end

@variables x S_0 [latexwrapper = identity] r d σ T δ Z W_τ [latexwrapper = identity]
const mu_lrm = log(S_0) + (r - d - σ^2 / 2) * δ^2
const s_lrm = σ * δ

function lrm_scores()
    lines = ["\\partial_{$name} \\log f = " * latex_weight(mu_lrm, s_lrm, tgt_var, 1, x, Z) for (name, tgt_var) in (("S_0", S_0), ("r", r), ("d", d), ("\\sigma", σ))]
    # T only enters through δ = sqrt(T), and dδ/dT = 1 / (2δ).
    score_T = simplify_fractions(generate_weight(mu_lrm, s_lrm, δ, 1, x, Z) / (2 * δ))
    return display_math([lines; "\\partial_{T} \\log f = " * latex_weight(score_T)])
end
lrm_delta() = display_math(["\\Delta = \\partial_{S_0} V_0 = \\mathbb{E}\\left(\\Phi(X_T) " * latex_weight(mu_lrm, s_lrm, S_0, 1, x, Z) * "\\right)"])
lrm_delta_wt() = display_math([replace(generate_latex_expression_base(S_0, r, d, σ, T, S_0, 1), "V \\partial" => "\\partial V")])
sixth_order() = display_math([generate_latex_expression_base(S_0, r, d, σ, T, r, 6)])

# Black and Scholes over the remaining time h: the derivatives in S of the call and of the binary call.
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
delta_call() = generate_delta(S * exp(-d * h) * N(D1) - K * exp(-r * h) * N(D2))
delta_binary() = substitute(expand_derivatives(Differential(S)(exp(-r * h) * N(D2))), Dict(Differential(S)(D2) => 1 / (S * σ * sqrt(h))))
conditional_derivatives() = display_math(["\\partial_{S} C_{\\text{call}} = " * latex_delta(delta_call()), "\\partial_{S} C_{\\text{binary}} = " * latex_delta(delta_binary())])
@variables λ κ
with_jumps(delta) = substitute(delta, Dict(d => d + λ * κ))
conditional_saltando_derivatives() = display_math(["\\partial_{S} C_{\\lambda \\kappa,\\text{call}} = " * latex_delta(with_jumps(delta_call())), "\\partial_{S} C_{\\lambda \\kappa,\\text{binary}} = " * latex_delta(with_jumps(delta_binary()))])

@variables μ s
gaussian_scores() = display_math(["\\partial_{\\mu} \\log f = " * latex_weight(μ, s, μ, 1, x, Z), "\\partial_{s} \\log f = " * latex_weight(μ, s, s, 1, x, Z)])

const mu_bs = log(S_0) + (r - d - σ^2 / 2) * T + σ * W_τ
const s_bs = σ * δ
weights(mu) = display_math(["w_{\\Delta} = " * latex_weight(mu, s_bs, S_0, 1, x, Z), "w_{\\Gamma} = " * latex_weight(mu, s_bs, S_0, 2, x, Z), "w_{\\mathcal{V}} = " * latex_weight(mu, s_bs, σ, 1, x, Z)])
vibrato_weights() = weights(mu_bs)
@variables J
saltando_weights() = weights(log(S_0) + (r - d - σ^2 / 2 - λ * κ) * T + σ * W_τ + J)
@variables β ψ
function subordinated_weights()
    mu_sub = log(S_0) + (r - d - ψ) * T + β * δ^2
    return display_math(["w_{\\Delta} = " * latex_weight(mu_sub, s_bs, S_0, 1, x, Z), "w_{\\Gamma} = " * latex_weight(mu_sub, s_bs, S_0, 2, x, Z)])
end
