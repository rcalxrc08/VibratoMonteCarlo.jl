# This file was generated, do not modify it. # hide
#hideall
@variables J λ κ
mu_saltando = log(S_0) + (r - d - σ^2 / 2 - λ * κ) * T + σ * W_τ + J
println("\$\$w_{\\Delta} = ", latex_weight(mu_saltando, s_bs, S_0, 1, x, Z), "\$\$")
println("\$\$w_{\\Gamma} = ", latex_weight(mu_saltando, s_bs, S_0, 2, x, Z), "\$\$")
println("\$\$w_{\\mathcal{V}} = ", latex_weight(mu_saltando, s_bs, σ, 1, x, Z), "\$\$")