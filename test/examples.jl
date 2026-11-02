@isdefined(TEST_UTILS_LOADED) || include("utils.jl")
include(joinpath(@__DIR__, "..", "etc", "readme.jl"))

@testset "examples" begin
    # README.md shows the code of the examples: after changing one, run julia etc/readme.jl.
    @testset "README is up to date" begin
        readme = replace(read(README_PATH, String), "\r\n" => "\n")
        @test readme_with_examples(readme) == readme
    end

    @testset "hello_world.jl" begin
        example = Module()
        redirect_stdout(devnull) do
            Base.include(example, joinpath(@__DIR__, "..", "examples", "hello_world.jl"))
        end
        # Pure LRM with a fine grid: relative errors observed 1.2e-7 (value), 4.5e-9 (delta).
        ref = blsprice(dual(100.0, 1.0), 100.0, 0.02, 1.0, 0.2, 0.01)
        @test DualNumbers.value(example.EuPrice) ≈ DualNumbers.value(ref) rtol = 1e-6
        @test DualNumbers.epsilon(example.EuPrice) ≈ DualNumbers.epsilon(ref) rtol = 1e-6
    end
end
