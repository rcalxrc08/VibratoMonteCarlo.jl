module VibratoMonteCarloTaylorSeriesExt
using VibratoMonteCarlo, TaylorSeries
import VibratoMonteCarlo: v_value

v_value(x::Taylor1) = v_value(x[0])
v_value(x::AbstractSeries) = @views x[0][1]

end # module
