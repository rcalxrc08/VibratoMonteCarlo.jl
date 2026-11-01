# This file was generated, do not modify it. # hide
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