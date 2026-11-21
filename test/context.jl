# The layered context (CLAUDE.md, CONTEXT.md, the stage contracts, references/) must follow the code: every file a
# contract lists exists, every file the module includes and every reference is reached from a contract, and every
# contract is routed from the root CONTEXT.md (a reference may instead be named by CLAUDE.md, as a rule for every change),
# and state/STATUS.md (layer 4) has its sections and an Updated date.
# A moved or added file without its contract fails here instead of leaving the map stale. Untracked files under src/ are
# not checked: only the files the module includes.
using Test

@testset "Context layers" begin
    root = normpath(joinpath(@__DIR__, ".."))
    rel(path) = replace(relpath(path, root), '\\' => '/')
    contracts = vcat([joinpath(root, d, "CONTEXT.md") for d in ("ext", "test", "docs", "bench", "etc")], [joinpath(root, "src", d, "CONTEXT.md") for d in readdir(joinpath(root, "src")) if isdir(joinpath(root, "src", d)) && occursin(r"^\d\d_", d)])
    routing = read(joinpath(root, "CONTEXT.md"), String)
    layer0 = read(joinpath(root, "CLAUDE.md"), String)
    listed = Set{String}()
    for contract in contracts
        @test isfile(contract)
        @test occursin("`$(rel(contract))`", routing)
        dir = dirname(contract)
        for line in eachline(contract)
            startswith(line, "| `") || continue
            first_cell = split(line, '|')[2]
            for m in eachmatch(r"`([^`*]+\.(?:jl|md))`", first_cell)
                name = m.captures[1]
                path = isfile(joinpath(dir, name)) ? joinpath(dir, name) : joinpath(root, name)
                @test isfile(path) || error("$(rel(contract)) lists $name, which does not exist")
                push!(listed, rel(path))
            end
        end
    end
    module_file = joinpath(root, "src", "VibratoMonteCarlo.jl")
    for m in eachmatch(r"include\(\"([^\"]+)\"\)", read(module_file, String))
        path = "src/" * m.captures[1]
        @test path in listed || error("$path is in no stage contract")
    end
    for file in readdir(joinpath(root, "references"))
        path = "references/$file"
        @test path in listed || occursin(path, layer0) || error("$path is reached from no contract nor CLAUDE.md")
    end
    status = joinpath(root, "state", "STATUS.md")
    @test isfile(status) || error("state/STATUS.md is missing")
    text = read(status, String)
    for section in ("## Current focus", "## Active work", "## Known issues")
        @test occursin(section, text) || error("state/STATUS.md has no '$section' section")
    end
    @test occursin(r"Updated: \d{4}-\d{2}-\d{2}", text) || error("state/STATUS.md has no 'Updated: YYYY-MM-DD' date")
end
