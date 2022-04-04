using .ForwardDiff

v_value(x::ForwardDiff.Dual{Tg,T,N}) where {Tg,T<:Real,N} = ForwardDiff.value(x)
v_mod(x::ForwardDiff.Dual{Tg,T,N}) where {Tg,T<:Real,N} = collect(ForwardDiff.partials(x))