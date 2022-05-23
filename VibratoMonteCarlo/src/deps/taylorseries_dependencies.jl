using .TaylorSeries

import Base.floor;
floor(x::Taylor1) =@views floor(x[0])
floor(x::AbstractSeries) =@views floor(x[0][1])
v_value(x::Taylor1) =@views x[0]
v_value(x::AbstractSeries) =@views x[0][1]
