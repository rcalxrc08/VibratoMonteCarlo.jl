using .TaylorSeries

v_value(x::Taylor1) =@views x[0]
v_value(x::AbstractSeries) =@views x[0][1]

# function BSplineKit.interpolate(
        # x::AbstractVector, y::AbstractVector{T1}, k::BSplineOrder,
        # ::Nothing = nothing,
    # ) where {T1 <: TaylorSeries.AbstractSeries{Float64}}
    # t = BSplineKit.SplineInterpolations.make_knots(x, order(k))
    # B = BSplineBasis(k, t; augment = Val(false))  # it's already augmented!

    # # If input data is integer, convert the spline element type to float.
    # # This also does the right thing when eltype(y) <: StaticArray.
    # T = eltype(y)

    # itp = SplineInterpolation(undef, B, x, T)
    # interpolate!(itp, y)
# end