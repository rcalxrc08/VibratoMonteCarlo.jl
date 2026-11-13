# Figures of the documentation, computed with the package at build time and saved as SVG in src/assets/figures, in a
# light and a dark variant: the stylesheet assets/figures.css shows the one of the theme of the reader. The theme is the
# one of the FinancialFFT.jl documentation.
using VibratoMonteCarlo, FinancialMonteCarlo, FinancialFFT, FinancialToolbox, ForwardDiff, HyperDualNumbers, Plots
# GR draws without a display on the CI runner
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
gr()

const FIGURES = mkpath(joinpath(@__DIR__, "src", "assets", "figures"))

# Categorical slots in fixed order (validated for colour vision deficiencies on adjacent pairs, in both modes), and
# the ink of the text, of the axes and of the grid.
const THEMES = (light = (series = ["#2563d9", "#0a8a6f", "#b86e0b", "#0b7fb3", "#2f8a3c", "#c2456e"], text = "#1f2328", muted = "#5b6470", axis = "#b9b4ab", grid = "#e6e3de", diverging = ["#2563d9", "#f7f6f3", "#c47a12"]), dark = (series = ["#4f8cff", "#00d4aa", "#ffb84d", "#38cfff", "#7ee787", "#ff6b9d"], text = "#f5f7fa", muted = "#a0aec0", axis = "#3a4466", grid = "#252d47", diverging = ["#4f8cff", "#151b2e", "#ffb84d"]))
# The secondary encoding of the series: the light slots 3-5 are below 3:1 contrast on the surface.
const MARKERS = [:circle, :rect, :utriangle, :diamond, :dtriangle, :star5]

function base_plot(theme; kwargs...)
    return plot(; size = (760, 430), background_color = :transparent, foreground_color_text = theme.muted, foreground_color_guide = theme.muted, foreground_color_axis = theme.axis, foreground_color_border = theme.axis, foreground_color_legend = nothing, legend_background_color = :transparent, legendfontcolor = theme.text, titlefontcolor = theme.text, gridcolor = theme.grid, gridalpha = 1, gridlinewidth = 1, framestyle = :axes, fontfamily = "sans-serif", titlefontsize = 12, guidefontsize = 10, tickfontsize = 9, legendfontsize = 9, titlelocation = :left, left_margin = 4Plots.mm, bottom_margin = 4Plots.mm, kwargs...)
end

function save_themes(draw, name)
    for (theme_name, theme) in pairs(THEMES)
        savefig(draw(theme), joinpath(FIGURES, "$name-$theme_name.svg"))
    end
end

const r, d, T, K, σ = 0.02, 0.0, 1.0, 100.0, 0.2
const zero_rate = ZeroRate(r)
value_delta_gamma(x::Hyper) = (x.value, x.epsilon1, x.epsilon12)
hyper_of(x) = hyper(x, 1.0, 1.0, 0.0)

# Delta and gamma of a binary call, the payoff with a jump at the strike: plain Monte Carlo has zero pathwise derivative,
# Vibrato gets them from hyper dual numbers, on 10^4 paths, against the formula of Black and Scholes.
function binary_greeks_figure(theme)
    spots = collect(range(70.0, 130.0, length = 25))
    fine = range(70.0, 130.0, length = 301)
    mc = MonteCarloConfiguration(10_000, 10; sampling = FinancialMonteCarlo.AntitheticMC())
    option = BinaryEuropeanOption(T, K)
    exact(S0) = value_delta_gamma(blsbin(hyper_of(S0), K, r, T, σ, d))
    vibrato(S0) = value_delta_gamma(pricer(BlackScholesProcess(σ, Underlying(hyper_of(S0), d)), zero_rate, Vibrato(mc, VibratoMonteCarloGaussLegendre()), option))
    estimates = vibrato.(spots)
    panels = map(((2, "delta"), (3, "gamma"))) do (k, name)
        p = base_plot(theme; title = "$name of a binary call, K = 100, T = 1", xlabel = "S₀", ylabel = name, legend = k == 2 ? :topright : :bottomleft)
        plot!(p, fine, [exact(S0)[k] for S0 in fine]; label = "Black and Scholes", color = theme.muted, linewidth = 1.5, linestyle = :dash)
        scatter!(p, spots, [e[k] for e in estimates]; label = "Vibrato, 10⁴ paths", color = theme.series[1], marker = MARKERS[1], markersize = 5, markerstrokewidth = 0)
        return p
    end
    return plot(panels...; layout = (1, 2), size = (1000, 400))
end

# The inner integration alone: PureLRM on Black and Scholes is one integral against the density at maturity, so its error
# is the error of the quadrature, here on the gamma of a binary call, whose payoff jumps at the strike.
function quadrature_figure(theme)
    option = BinaryEuropeanOption(T, K)
    model = BlackScholesProcess(σ, Underlying(hyper_of(K), d))
    mc = MonteCarloConfiguration(1, 1)
    exact = value_delta_gamma(blsbin(hyper_of(K), K, r, T, σ, d))[3]
    floor_error = 1e-16
    error(vb) = max(abs(value_delta_gamma(pricer(model, zero_rate, PureLRM(mc, vb), option))[3] - exact) / abs(exact), floor_error)
    families = ["VibratoMonteCarloGaussLegendre(n), split at the strike" => [(2n, VibratoMonteCarloGaussLegendre(n)) for n in [2, 4, 6, 8, 12, 16, 24, 32, 48]], "VibratoMonteCarlo.VibratoMonteCarloGaussHermite(n)" => [(n, VibratoMonteCarlo.VibratoMonteCarloGaussHermite(n)) for n in [4, 8, 16, 32, 64, 128, 200]], "VibratoMonteCarloAnalytic(n), uniform grid on ±7" => [(n, VibratoMonteCarloAnalytic(n)) for n in [8, 16, 32, 64, 128, 256, 512, 1024]], "VibratoMonteCarloStandard(n), Monte Carlo" => [(n, VibratoMonteCarloStandard(n)) for n in [8, 32, 128, 512, 2048, 8192]]]
    p = base_plot(theme; title = "Relative error of the gamma of a binary call (K = S₀ = 100, T = 1), one inner integral", xlabel = "integration nodes", ylabel = "relative error", xscale = :log10, yscale = :log10, legend = :outerright, size = (1000, 430), ylims = (1e-16, 10), yticks = 10.0 .^ (-16:2:0))
    for (j, (name, family)) in enumerate(families)
        plot!(p, first.(family), [error(vb) for (_, vb) in family]; label = name, color = theme.series[j], linewidth = 2, marker = MARKERS[j], markersize = 5, markerstrokewidth = 0)
    end
    return p
end

# Variance Gamma: the FFT density of the last step is singular when dt / κ is small, while VibratoSaltando conditions on
# the subordinator and does not depend on the time grid. Error of the gamma of a call against the Lewis formula.
function variance_gamma_figure(theme)
    κ = 0.2
    model = VarianceGammaProcess(0.2, -0.15, κ, Underlying(hyper_of(100.0), d))
    option = EuropeanOption(T, K)
    exact = value_delta_gamma(pricer(model, zero_rate, LewisMethod(400.0, 200_000), option))[3]
    steps = [2, 4, 8, 16, 32]
    mc(n) = MonteCarloConfiguration(20_000, n; sampling = FinancialMonteCarlo.AntitheticMC())
    rel(x) = abs(value_delta_gamma(x)[3] - exact) / abs(exact)
    fft = [rel(pricer(model, zero_rate, Vibrato(mc(n), VibratoMonteCarloAnalytic()), option)) for n in steps]
    saltando = [rel(pricer(model, zero_rate, VibratoSaltando(mc(n), VibratoMonteCarloGaussLegendre()), option)) for n in steps]
    ratio = (T ./ steps) ./ κ
    p = base_plot(theme; title = "Variance Gamma call (κ = 0.2, T = 1): relative error of the gamma, 2 × 10⁴ paths", xlabel = "dt / κ of the last step", ylabel = "relative error", xscale = :log10, yscale = :log10, legend = :outerright, size = (1000, 430), xflip = true, xticks = (ratio, string.(round.(ratio, sigdigits = 3))), ylims = (1e-4, 1e5), yticks = 10.0 .^ (-4:2:4))
    plot!(p, ratio, fft; label = "Vibrato, FFT density of the last step", color = theme.series[1], linewidth = 2, marker = MARKERS[1], markersize = 5, markerstrokewidth = 0)
    plot!(p, ratio, saltando; label = "VibratoSaltando, conditioned on G_T", color = theme.series[2], linewidth = 2, marker = MARKERS[2], markersize = 5, markerstrokewidth = 0)
    return p
end

save_themes(binary_greeks_figure, "binary_greeks")
save_themes(quadrature_figure, "quadrature")
save_themes(variance_gamma_figure, "variance_gamma")
# The hero of the home page: the gamma of a binary call across spot and maturity, from PureLRM with Gauss-Legendre (one
# integral against the density at maturity, exact up to rounding) differentiated with hyper-dual numbers, and Vibrato
# Monte Carlo estimates of the same gamma as points; drawn in perspective by surface.jl.
include(joinpath(@__DIR__, "surface.jl"))
let spots = collect(range(60.0, 140.0; length = 46)), maturities = collect(range(0.15, 2.0; length = 34))
    lrm = PureLRM(MonteCarloConfiguration(1, 1), VibratoMonteCarloGaussLegendre())
    gamma(method, S, t) = pricer(BlackScholesProcess(σ, Underlying(hyper(S, 1.0, 1.0, 0.0), d)), zero_rate, method, BinaryEuropeanOption(t, K)).epsilon12
    gammas = [gamma(lrm, S, t) for t in maturities, S in spots]
    vibrato = Vibrato(MonteCarloConfiguration(100_000, 10; sampling = FinancialMonteCarlo.AntitheticMC()), VibratoMonteCarloGaussLegendre())
    points = [(S, t, gamma(vibrato, S, t), true) for t in (0.5, 1.5) for S = 70.0:10.0:130.0]
    m = ceil(maximum(abs, gammas) / 0.001) * 0.001
    write_surface(FIGURES, "hero-surface", spots, maturities, gammas; zlims = (-m, m), points = points, xticks = [s => string(round(Int, s)) for s = 60:20:140], yticks = [0.5 => "6M", 1.0 => "1Y", 1.5 => "18M", 2.0 => "2Y"], zticks = [g => string(g) for g = (-m):0.002:m], xlabel = "spot", ylabel = "maturity", zlabel = "gamma", label = "Gamma of a binary call across spot and maturity")
end
