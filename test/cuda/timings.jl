# Timings of the parts moved to the GPU (#51), on the device and, where it helps, on the host. Not a test: run with
#     julia --project=test/cuda test/cuda/timings.jl
# BENCH_PART=sub|jumps|fft|am runs one group. Each line is the best of 3 runs after a warm-up, synchronized with the device.
using CUDA, VibratoMonteCarlo, FinancialMonteCarlo, FinancialFFT, HyperDualNumbers
const VMC = VibratoMonteCarlo
CUDA.allowscalar(false)
h(x) = hyper(x, 1.0, 1.0, 0.0)
rf = ZeroRate(0.02);
eu = EuropeanOption(1.0, 100.0)
gl = VibratoMonteCarloGaussLegendre()
gpu(s = FinancialMonteCarlo.StandardMC(); n = 100_000, steps = 10) = MonteCarloConfiguration(n, steps; sampling = s, mode = FinancialMonteCarlo.CudaMode())
cpu(s = FinancialMonteCarlo.StandardMC(); n = 100_000, steps = 10) = MonteCarloConfiguration(n, steps; sampling = s)
function t(label, f)
    f()
    CUDA.synchronize()
    best = minimum(1:3) do _
        CUDA.@elapsed begin
            f()
        end
    end
    println(rpad(label, 46), round(best, digits = 4), " s")
end
vg = VarianceGammaProcess(0.2, -0.15, 0.2, Underlying(h(100.0), 0.0))
nig = NormalInverseGaussianProcess(0.2, -0.1, 0.2, Underlying(h(100.0), 0.0))
kou = KouProcess(0.2, 5.0, 0.3, 30.0, 20.0, Underlying(h(100.0), 0.0))
mer = MertonProcess(0.2, 5.0, -0.05, 0.1, Underlying(h(100.0), 0.0))
part = get(ENV, "BENCH_PART", "all")
if part in ("all", "sub")
    for (n, m) in (("VG", vg), ("NIG", nig)), (sn, s) in (("plain", FinancialMonteCarlo.StandardMC()),)
        t("G_T $n $sn 1M (simulate_subordinator)", () -> VMC.simulate_subordinator(m, gpu(s; n = 1_000_000), 1.0))
        t("VibratoSaltando $n $sn 100k GL (GPU)", () -> pricer(m, rf, VibratoSaltando(gpu(s), gl), eu))
    end
end
if part in ("all", "jumps")
    for (n, m) in (("Kou", kou), ("Merton", mer)), (sn, s) in (("plain", FinancialMonteCarlo.StandardMC()), ("anti", FinancialMonteCarlo.AntitheticMC()))
        t("jump times $n $sn 100k×10 (simulate_log…)", () -> VMC.simulate_log_with_jump_times(m, rf, gpu(s), 1.0))
        t("VibratoSaltando $n $sn 100k GL (GPU)", () -> pricer(m, rf, VibratoSaltando(gpu(s), gl), eu))
    end
end
if part in ("all", "fft")
    an = VibratoMonteCarloAnalytic()
    for (n, m) in (("Kou", kou), ("NIG", nig))
        t("FFT Vibrato $n 20k×10 (GPU)", () -> pricer(m, rf, Vibrato(gpu(; n = 20_000), an), eu))
        t("FFT Vibrato $n 20k×10 (CPU)", () -> pricer(m, rf, Vibrato(cpu(; n = 20_000), an), eu))
        t("PureLRM $n (GPU)", () -> pricer(m, rf, PureLRM(gpu(), an), eu))
    end
end
if part in ("all", "am")
    bs = BlackScholesProcess(0.2, Underlying(h(100.0), 0.0))
    hes = HestonProcess(0.2, 0.2, 0.0, 0.2, 0.2, 0.2, Underlying(h(100.0), 0.0))
    put = AmericanOption(1.0, 100.0, false)
    an = VibratoMonteCarloAnalytic()
    t("American BS put 20k×10 GL (GPU)", () -> pricer(bs, rf, Vibrato(gpu(; n = 20_000), gl), put))
    t("American Heston put 20k×10 GL (GPU)", () -> pricer(hes, rf, Vibrato(gpu(; n = 20_000), gl), put))
    t("American Kou put 20k×10 Analytic (GPU)", () -> pricer(kou, rf, Vibrato(gpu(; n = 20_000), an), put))
    t("American Kou put 20k×10 Analytic (CPU)", () -> pricer(kou, rf, Vibrato(cpu(; n = 20_000), an), put))
    t("PureLRM BS GL (GPU)", () -> pricer(bs, rf, PureLRM(gpu(), gl), eu))
end
