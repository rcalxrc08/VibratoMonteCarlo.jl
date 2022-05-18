using .TaylorSeries

import Base.floor;
floor(x::AbstractSeries) = floor(x[0])
v_value(x::AbstractSeries) = x[0]

function bell_n(prev_bell, X)
    if length(X) == 0
        return 1.0
    end
    if length(X) == 1
        return X[1]
    end
    if length(X) == 2
        return X[1]^2 + X[2]
    end
    result = 0.0
    n = length(X) - 1
    # for i in 0:n
    # coeff_bin_i=binomial(n,i);
    # @views X_i_1=X[i+1]
    # @views result+=coeff_bin_i*prev_bell[n-i+1]*X_i_1
    # end
    return @views sum(binomial(n, i) * prev_bell[n-i+1] * X[i+1] for i = 0:n)
end

function faa_di_bruno_bell(x::Taylor1)
    max_ord = get_order(x)
    orders = 1:max_ord
    @views derivatives = [differentiate(x, i)[0] for i in orders]
    prev_bell = Array{Float64}(undef, max_ord + 1)
    @views prev_bell[1] = 1.0
    for i in orders
        @views result = bell_n(prev_bell, derivatives[1:i])
        @views prev_bell[i+1] = result
    end
    @. prev_bell = prev_bell / factorial(0:max_ord)
    return Taylor1(prev_bell)
end

function v_mod(x::AbstractSeries)
    return faa_di_bruno_bell(x)
end