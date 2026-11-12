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
merton(S) = MertonProcess(0.2, 5.0, -0.05, 0.1, Underlying(S, 0.0))
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

    # Plain MC draws G_T on the device (FinancialMonteCarlo's Gamma and inverse Gaussian kernel), so only its law can be
    # compared with the host: the empirical distribution function at the deciles, within 4 standard errors.
    @testset "G_T drawn on the device" begin
        n = 200_000
        for model in (var_gamma(100.0), nig(100.0))
            law = FinancialMonteCarlo.compute_subordinator_law(model, model.κ, 1.0)
            (G_T, weight) = VMC.simulate_subordinator(model, gpu(; n = n), 1.0)
            @test G_T isa CuArray
            @test weight === true
            G = Array(G_T)
            for p = 0.1:0.1:0.9
                @test abs(count(<=(VMC.Distributions.quantile(law, p)), G) / n - p) < 4 * sqrt(p * (1 - p) / n)
            end
            # set_seed! makes the draws reproducible (the kernel is seeded from the host generator), a new seed changes them
            draw(mc) = (FinancialMonteCarlo.set_seed!(mc); Array(VMC.simulate_subordinator(model, mc, 1.0)[1]))
            @test draw(gpu(; n = n)) == draw(gpu(; n = n))
            @test draw(gpu(; n = n)) != draw(MonteCarloConfiguration(n, 10; mode = FinancialMonteCarlo.CudaMode(), seed = 7))
        end
        # the weight of κ and T on the device, against the host on the same G_T: ≤ 1e-15 observed
        model = var_gamma(100.0; κ = hyper_of(0.2))
        (G_T, weight) = VMC.simulate_subordinator(model, gpu(; n = 10_000), hyper_of(1.0))
        @test weight isa CuArray
        host_weight = map(g -> VMC.v_mod(VMC.subordinator_log_density(model, g, hyper_of(1.0), model.κ)), Array(G_T))
        @test all(isapprox.(value_and_derivatives.(Array(weight)), value_and_derivatives.(host_weight); rtol = 1e-12))
    end

    # AntitheticMC draws the quantiles of G_T on the host and uploads them, so the device gives the host result: ≤ 1.5e-15
    # observed. The antithetic halves are ordered as FinancialMonteCarlo's path_sides. (FinancialMonteCarlo refuses SobolMC
    # in CudaMode.)
    @testset "antithetic G_T as the host" begin
        anti = FinancialMonteCarlo.AntitheticMC()
        cases = ((var_gamma(hyper_of(100.0)), eu), (var_gamma(100.0; κ = hyper_of(0.2)), eu), (nig(hyper_of(100.0)), bin))
        for (model, opt) in cases
            @test value_and_derivatives(pricer(model, rf, VibratoSaltando(gpu(anti), gl), opt)) ≈ value_and_derivatives(pricer(model, rf, VibratoSaltando(cpu(anti), gl), opt)) rtol = 1e-12
            @test value_and_derivatives(pricer(model, rf, VMC.ConditionalSaltando(gpu(anti)), opt)) ≈ value_and_derivatives(pricer(model, rf, VMC.ConditionalSaltando(cpu(anti)), opt)) rtol = 1e-12
        end
        for model in (var_gamma(100.0), nig(100.0))
            law = FinancialMonteCarlo.compute_subordinator_law(model, model.κ, 1.0)
            mc = gpu(anti; n = 10_000)
            G = Array(VMC.simulate_subordinator(model, mc, 1.0)[1])
            sides = FinancialMonteCarlo.path_sides(reshape(1:mc.Nsim, :, 1), mc)
            first_half, second_half = (vec(collect(side[1])) for side in sides)
            @test VMC.Distributions.cdf.(law, G[first_half]) .+ VMC.Distributions.cdf.(law, G[second_half]) ≈ ones(length(first_half)) rtol = 1e-8
        end
        device = pricer(nig(hyper_of(100.0)), rf, VibratoSaltando(gpu(anti), gl), smile)
        host = pricer(nig(hyper_of(100.0)), rf, VibratoSaltando(cpu(anti), gl), smile)
        @test reduce(vcat, value_and_derivatives.(device)) ≈ reduce(vcat, value_and_derivatives.(host)) rtol = 1e-12
        # a smile with its strikes on the device, in either mode: its options are built on the host
        cu_smile = FinancialMonteCarlo.EuropeanOptionSmile(1.0, CuArray([90.0, 100.0, 115.0]), CuArray([true, false, true]))
        @test reduce(vcat, value_and_derivatives.(pricer(nig(hyper_of(100.0)), rf, VibratoSaltando(gpu(anti), gl), cu_smile))) ≈ reduce(vcat, value_and_derivatives.(device)) rtol = 1e-12
        @test reduce(vcat, value_and_derivatives.(pricer(nig(hyper_of(100.0)), rf, VibratoSaltando(cpu(anti), gl), cu_smile))) ≈ reduce(vcat, value_and_derivatives.(host)) rtol = 1e-12
    end

    # The last jumps gathered on the device, against the old host algorithm on the same draws: the CURAND diffusion again,
    # taken to the host, the same JumpList of FinancialMonteCarlo (same host generator state), and path by path its jumps
    # in time order added to the row, the grid point before the step of the last jump plus the jumps of that step.
    @testset "last jumps on the device" begin
        for model in (kou(hyper_of(100.0)), merton(hyper_of(100.0))), sampling in (FinancialMonteCarlo.StandardMC(), FinancialMonteCarlo.AntitheticMC())
            mc = gpu(sampling; n = 20_000)
            FinancialMonteCarlo.set_seed!(mc)
            (T_jumps, X_f) = VMC.simulate_log_with_jump_times(model, rf, mc, 1.0)
            @test T_jumps isa CuArray && X_f isa CuArray
            FinancialMonteCarlo.set_seed!(mc)
            drift = rf.r - VMC.drift_correction(model)
            X = FinancialMonteCarlo.get_matrix_type(mc, model, typeof(drift))
            FinancialMonteCarlo.simulate!(X, BrownianMotion(model.σ, drift), mc, 1.0)
            X = Array(X)
            jumps = FinancialMonteCarlo.JumpList{Float64}()
            FinancialMonteCarlo.add_jumps!(jumps, model, mc, 1.0)
            offsets, steps, sizes = FinancialMonteCarlo.jumps_by_path(jumps, mc.Nsim)
            T_ref = zeros(mc.Nsim)
            X_ref = zeros(mc.Nsim)
            for i = 1:mc.Nsim
                row = X[i, :]
                path_jumps = offsets[i]:(offsets[i+1]-1)
                last_idx = 0
                jumps_in_step = 0.0
                for k in path_jumps[sortperm(steps[path_jumps])]
                    row[steps[k]:end] .+= sizes[k]
                    jumps_in_step = steps[k] == last_idx ? jumps_in_step + sizes[k] : sizes[k]
                    last_idx = steps[k]
                    T_ref[i] = (last_idx - 2) / mc.Nstep
                    X_ref[i] = row[last_idx-1] + jumps_in_step
                end
            end
            @test count(>(0), T_ref) > mc.Nsim / 2
            @test Array(T_jumps) ≈ T_ref rtol = 1e-14
            @test all(all(isapprox.(value_and_derivatives(x), value_and_derivatives(log(model.underlying.S0) + y); rtol = 1e-12)) for (x, y) in zip(Array(X_f), X_ref))
        end
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
        # G_T from the device sampler: worst relative errors over 6 seeds (value, delta, gamma) 1.2e-3, 1.7e-3, 1.5e-3.
        for model in (var_gamma(hyper_of(100.0)), nig(hyper_of(100.0))), method in (mc -> VibratoSaltando(mc, gl), VMC.ConditionalSaltando)
            ref = value_and_derivatives(pricer(model, rf, lewis, eu))
            @test all(abs.(value_and_derivatives(pricer(model, rf, method(gpu()), eu)) .- ref) .< [2.5e-3, 3.5e-3, 3e-3] .* abs.(ref))
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

    # The FFT density of the last step (Vibrato and PureLRM on non-Ito processes) is computed and evaluated on the device.
    @testset "FFT density on the device" begin
        an = VibratoMonteCarloAnalytic()
        for model in (kou(hyper_of(100.0)), nig(hyper_of(100.0)), var_gamma(hyper_of(100.0); κ = 0.1))
            # the FFTs of the device and of the host: values within 2.1e-15, the moments of the grid within 1e-14
            host_law = VMC.fft_step_law(model, rf, 0.25)
            device_law = VMC.fft_step_law(model, rf, 0.25, FinancialMonteCarlo.CudaMode())
            @test device_law.density.values isa CuArray
            @test maximum(abs.(Array(device_law.density.values) .- host_law.density.values)) < 1e-13
            @test device_law.mean ≈ host_law.mean atol = 1e-13
            @test device_law.std ≈ host_law.std rtol = 1e-12
            # the integration of the same paths with the same law (the device one, also taken to the host), device against
            # host: values ≤ 1.2e-13, gammas ≤ 1.1e-12 apart (absolute)
            same_law = VMC.FFTStepLaw(FinancialFFT.adapt(Array, device_law.density), device_law.mean, device_law.std)
            S = hyper.(100.0 .* exp.(0.2 .* randn(3_000)), 1.0, 1.0, 0.0)
            for opt in (eu, bin, smile)
                host = VMC.aug_paths(same_law, opt, an, S, cpu())
                host = host isa AbstractVector{<:AbstractVector} ? permutedims(reduce(hcat, host)) : host
                device = Array(VMC.aug_paths(device_law, opt, an, CuArray(S), gpu()))
                @test all(all(isapprox.(value_and_derivatives(x), value_and_derivatives(y); rtol = 1e-12, atol = 1e-11)) for (x, y) in zip(device, host))
            end
            # PureLRM has no paths: the device gives the host result up to the FFTs
            for opt in (eu, bin)
                @test value_and_derivatives(pricer(model, rf, PureLRM(gpu(), an), opt)) ≈ value_and_derivatives(pricer(model, rf, PureLRM(cpu(), an), opt)) rtol = 1e-9
            end
            @test reduce(vcat, value_and_derivatives.(pricer(model, rf, PureLRM(gpu(), an), smile))) ≈ reduce(vcat, value_and_derivatives.(pricer(model, rf, PureLRM(cpu(), an), smile))) rtol = 1e-9
        end
        # parameters with derivatives in the density computed on the device, and nested ForwardDiff duals there
        kou_σ(σ) = KouProcess(σ, 5.0, 0.3, 30.0, 20.0, Underlying(100.0, 0.0))
        @test value_and_derivatives(pricer(kou_σ(hyper_of(0.2)), rf, PureLRM(gpu(), an), eu)) ≈ value_and_derivatives(pricer(kou_σ(hyper_of(0.2)), rf, PureLRM(cpu(), an), eu)) rtol = 1e-9
        fft_vibrato(x) = pricer(KouProcess(x[2], 5.0, 0.3, 30.0, 20.0, Underlying(x[1], 0.0)), rf, Vibrato(gpu(; n = 5_000, steps = 4), an), eu)
        hess = ForwardDiff.hessian(fft_vibrato, [100.0, 0.2])
        @test hess[1, 1] ≈ fft_vibrato([hyper_of(100.0), 0.2]).epsilon12 rtol = 1e-9
        @test hess[2, 2] ≈ fft_vibrato([100.0, hyper_of(0.2)]).epsilon12 rtol = 1e-9
        @test hess[1, 2] ≈ fft_vibrato([hyper(100.0, 1.0, 0.0, 0.0), hyper(0.2, 0.0, 1.0, 0.0)]).epsilon12 rtol = 1e-9
        # CURAND paths, so against Lewis: worst relative errors over 6 seeds (value, delta, gamma)
        # 1.7e-2, 1.2e-2, 4.0e-3 on Kou, 1.3e-2, 9.0e-3, 2.5e-3 on NIG (50k paths, 4 steps)
        for model in (kou(hyper_of(100.0)), nig(hyper_of(100.0)))
            ref = value_and_derivatives(pricer(model, rf, lewis, eu))
            @test all(abs.(value_and_derivatives(pricer(model, rf, Vibrato(gpu(; n = 50_000, steps = 4), an), eu)) .- ref) .< [3.5e-2, 2.5e-2, 8e-3] .* abs.(ref))
        end
    end

    # The paths go to the host where the code indexes them (FFT density, exercise policy): they run, near the CPU ones.
    @testset "host parts" begin
        put = AmericanOption(1.0, 100.0, false)
        for model in (bs(hyper_of(100.0)), heston(hyper_of(100.0)))
            device = pricer(model, rf, Vibrato(gpu(; n = 20_000), gl), put)
            host = pricer(model, rf, Vibrato(cpu(; n = 20_000), gl), put)
            @test abs(device.value / host.value - 1) < 3e-2
        end
        @test pricer(bs(100.0), rf, PureLRM(gpu(), VibratoMonteCarloStandard(64)), eu) == pricer(bs(100.0), rf, PureLRM(cpu(), VibratoMonteCarloStandard(64)), eu)
    end
end
