using .TaylorSeries

v_value(x::Taylor1) =@views x[0]
v_value(x::AbstractSeries) =@views x[0][1]
