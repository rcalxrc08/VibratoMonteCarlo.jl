# This file was generated, do not modify it. # hide
#hideall
@variables x S_0 [latexwrapper = identity] r d σ T δ Z W_τ [latexwrapper = identity]
mu_bs = log(S_0) + (r - d - σ^2 / 2) * T + σ * W_τ
s_bs = σ * δ
println("\$\$w_{\\Delta} = ", latex_weight(mu_bs, s_bs, S_0, 1, x, Z), "\$\$")
println("\$\$w_{\\Gamma} = ", latex_weight(mu_bs, s_bs, S_0, 2, x, Z), "\$\$")
println("\$\$w_{\\mathcal{V}} = ", latex_weight(mu_bs, s_bs, σ, 1, x, Z), "\$\$")