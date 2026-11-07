# GPU tests, opt-in: they need an NVIDIA GPU and download the CUDA artifacts, so they are not part of runtests.jl. Run with
#     julia --project=test/cuda -e 'import Pkg; Pkg.instantiate(); include("test/cuda/runtests.jl")'
using CUDA, Test, VibratoMonteCarlo, FinancialMonteCarlo, FinancialFFT, FinancialToolbox, HyperDualNumbers, ForwardDiff
const VMC = VibratoMonteCarlo
CUDA.functional() || error("no functional CUDA device")
# Every GPU array is used whole: an indexing by scalars is a bug.
CUDA.allowscalar(false)

value_and_derivatives(x::Hyper) = [x.value, x.epsilon1, x.epsilon12]
hyper_of(x) = hyper(x, 1.0, 1.0, 0.0)
rf = ZeroRate(0.02)
eu = EuropeanOption(1.0, 100.0)
bin = BinaryEuropeanOption(1.0, 100.0)
smile = FinancialMonteCarlo.EuropeanOptionSmile(1.0, [90.0, 100.0, 115.0], [true, false, true])
bs(S, σ = 0.2) = BlackScholesProcess(σ, Underlying(S, 0.0))
kou(S) = KouProcess(0.2, 5.0, 0.3, 30.0, 20.0, Underlying(S, 0.0))
var_gamma(S; κ = 0.2) = VarianceGammaProcess(0.2, -0.15, κ, Underlying(S, 0.0))
nig(S) = NormalInverseGaussianProcess(0.2, -0.1, 0.2, Underlying(S, 0.0))
heston(S) = HestonProcess(0.2, 0.2, 0.0, 0.2, 0.2, 0.2, Underlying(S, 0.0))
gl = VibratoMonteCarloGaussLegendre()
lewis = LewisMethod(400.0, 200000)
gpu(method = FinancialMonteCarlo.StandardMC(); n = 100_000, steps = 10) = MonteCarloConfiguration(n, steps; sampling = method, mode = FinancialMonteCarlo.CudaMode())
cpu(method = FinancialMonteCarlo.StandardMC(); n = 100_000, steps = 10) = MonteCarloConfiguration(n, steps; sampling = method)

@testset "CUDA" begin
    # The matrix of paths × nodes on the device against the per path integration on the host: ≤ 3.9e-11 observed.
    @testset "inner integration, device against host" begin
        mu = hyper.(log(100.0) .+ 0.2 .* randn(2_000), 1.0, 1.0, 0.0)
        sigma = hyper(0.06, 0.0, 0.0, 0.0)
        for vb_mc in (gl, VibratoMonteCarloAnalytic(), VMC.VibratoMonteCarloGaussHermite(64)), opt in (eu, bin, smile)
            args = VMC.LRMInterfaceArguments(VMC.init_lrm_vec(vb_mc, cpu()), opt)
            host = VMC.lrm_paths(mu, sigma, args, cpu(), vb_mc)
            host = host isa AbstractVector{<:AbstractVector} ? permutedims(reduce(hcat, host)) : host
            device = Array(VMC.lrm_paths(CuArray(mu), sigma, args, gpu(), vb_mc))
            @test all(all(isapprox.(getfield.(device, f), getfield.(host, f); atol = 1e-9)) for f in (:value, :epsilon1, :epsilon12))
        end
    end

    # The same CURAND draws as FinancialMonteCarlo's GPU Heston, which floors the variance at 0 instead of eps: 1.4e-8 observed.
    @testset "Heston paths" begin
        mc = gpu(; n = 4_000)
        FinancialMonteCarlo.set_seed!(mc)
        (V, X) = VMC.simulate_with_vol(heston(100.0), rf, mc, 1.0)
        FinancialMonteCarlo.set_seed!(mc)
        @test Array(X) ≈ Array(FinancialMonteCarlo.simulate(heston(100.0), rf, mc, 1.0)) rtol = 1e-6
        @test V isa CuArray
    end

    # G_T is drawn on the host, so the device gives the host result: ≤ 1.5e-15 observed.
    @testset "same G_T as the host" begin
        cases = ((var_gamma(hyper_of(100.0)), eu), (var_gamma(100.0; κ = hyper_of(0.2)), eu), (nig(hyper_of(100.0)), bin))
        for (model, opt) in cases
            @test value_and_derivatives(pricer(model, rf, VibratoSaltando(gpu(), gl), opt)) ≈ value_and_derivatives(pricer(model, rf, VibratoSaltando(cpu(), gl), opt)) rtol = 1e-12
            @test value_and_derivatives(pricer(model, rf, VMC.ConditionalSaltando(gpu()), opt)) ≈ value_and_derivatives(pricer(model, rf, VMC.ConditionalSaltando(cpu()), opt)) rtol = 1e-12
        end
        device = pricer(nig(hyper_of(100.0)), rf, VibratoSaltando(gpu(), gl), smile)
        host = pricer(nig(hyper_of(100.0)), rf, VibratoSaltando(cpu(), gl), smile)
        @test reduce(vcat, value_and_derivatives.(device)) ≈ reduce(vcat, value_and_derivatives.(host)) rtol = 1e-12
        # a smile with its strikes on the device, in either mode: its options are built on the host
        cu_smile = FinancialMonteCarlo.EuropeanOptionSmile(1.0, CuArray([90.0, 100.0, 115.0]), CuArray([true, false, true]))
        @test reduce(vcat, value_and_derivatives.(pricer(nig(hyper_of(100.0)), rf, VibratoSaltando(gpu(), gl), cu_smile))) ≈ reduce(vcat, value_and_derivatives.(device)) rtol = 1e-12
        @test reduce(vcat, value_and_derivatives.(pricer(nig(hyper_of(100.0)), rf, VibratoSaltando(cpu(), gl), cu_smile))) ≈ reduce(vcat, value_and_derivatives.(host)) rtol = 1e-12
    end

    # CURAND paths, so only against the references. Worst relative errors over 6 GPU seeds with 100k paths (value, delta,
    # gamma): ≤ 1.2e-2, 5.3e-3, 1.0e-2 and 9.3e-2 for the gamma of the binary; the tolerances are about twice those.
    @testset "against references" begin
        anti = FinancialMonteCarlo.AntitheticMC()
        bs_ref = value_and_derivatives(blsprice(hyper_of(100.0), 100.0, 0.02, 1.0, 0.2, 0.0))
        tol = [2.5e-2, 1e-2, 2e-2]
        cases = (
            (pricer(bs(hyper_of(100.0)), rf, Vibrato(gpu(), gl), eu), bs_ref, tol),
            (pricer(bs(hyper_of(100.0)), rf, Vibrato(gpu(anti), gl), eu), bs_ref, tol),
            (pricer(bs(hyper_of(100.0)), rf, Vibrato(gpu(anti), VibratoMonteCarloStandard(64)), eu), bs_ref, tol),
            (pricer(bs(hyper_of(100.0)), rf, Vibrato(gpu(), gl), bin), value_and_derivatives(blsbin(hyper_of(100.0), 100.0, 0.02, 1.0, 0.2, 0.0)), [2.5e-2, 1e-2, 0.2]),
            (pricer(heston(hyper_of(100.0)), rf, Vibrato(gpu(), gl), eu), value_and_derivatives(pricer(heston(hyper_of(100.0)), rf, lewis, eu)), tol),
            (pricer(kou(hyper_of(100.0)), rf, VibratoSaltando(gpu(), gl), eu), value_and_derivatives(pricer(kou(hyper_of(100.0)), rf, lewis, eu)), tol),
            (pricer(kou(hyper_of(100.0)), rf, VMC.ConditionalSaltando(gpu()), eu), value_and_derivatives(pricer(kou(hyper_of(100.0)), rf, lewis, eu)), tol),
        )
        for (res, ref, tol) in cases
            @test all(abs.(value_and_derivatives(res) .- ref) .< tol .* abs.(ref))
        end
        smile_res = pricer(bs(hyper_of(100.0)), rf, Vibrato(gpu(), gl), smile)
        smile_ref = blsprice.(hyper_of(100.0), [90.0, 100.0, 115.0], 0.02, 1.0, 0.2, 0.0, [true, false, true])
        @test all(abs(r.value - q.value) < 2.5e-2 * q.value for (r, q) in zip(smile_res, smile_ref))
    end

    # Nested ForwardDiff duals on the device, on the same CURAND draws as hyper duals.
    @testset "ForwardDiff Hessian" begin
        hess = ForwardDiff.hessian(x -> pricer(bs(x[1], x[2]), rf, Vibrato(gpu(), gl), eu), [100.0, 0.2])
        @test hess[1, 1] ≈ pricer(bs(hyper_of(100.0)), rf, Vibrato(gpu(), gl), eu).epsilon12 rtol = 1e-10
        @test hess[2, 2] ≈ pricer(bs(100.0, hyper_of(0.2)), rf, Vibrato(gpu(), gl), eu).epsilon12 rtol = 1e-10
        @test hess[1, 2] ≈ pricer(bs(hyper(100.0, 1.0, 0.0, 0.0), hyper(0.2, 0.0, 1.0, 0.0)), rf, Vibrato(gpu(), gl), eu).epsilon12 rtol = 1e-10
    end

    # The paths go to the host where the code indexes them (FFT density, exercise policy): they run, near the CPU ones.
    @testset "host parts" begin
        fft = pricer(kou(hyper_of(100.0)), rf, Vibrato(gpu(; n = 2_000, steps = 2), VibratoMonteCarloAnalytic()), eu)
        @test abs(fft.value / pricer(kou(100.0), rf, lewis, eu) - 1) < 5e-2
        put = AmericanOption(1.0, 100.0, false)
        for model in (bs(hyper_of(100.0)), heston(hyper_of(100.0)))
            device = pricer(model, rf, Vibrato(gpu(; n = 20_000), gl), put)
            host = pricer(model, rf, Vibrato(cpu(; n = 20_000), gl), put)
            @test abs(device.value / host.value - 1) < 3e-2
        end
        @test pricer(bs(100.0), rf, PureLRM(gpu(), VibratoMonteCarloStandard(64)), eu) == pricer(bs(100.0), rf, PureLRM(cpu(), VibratoMonteCarloStandard(64)), eu)
    end
end
