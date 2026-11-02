# Copies the examples into README.md, between the markers
#     <!-- BEGIN examples/<name>.jl -->
#     <!-- END examples/<name>.jl -->
# GitHub cannot include a file in a README. Run from the repository root after changing an example:
#     julia etc/readme.jl
# test/examples.jl checks that README.md is up to date.

const README_PATH = joinpath(@__DIR__, "..", "README.md")

example_code(example) = rstrip(replace(read(joinpath(@__DIR__, "..", example), String), "\r\n" => "\n"))
readme_block(example) = "<!-- BEGIN $example -->\n```julia\n" * example_code(example) * "\n```\n<!-- END $example -->"

function readme_with_examples(readme)
    readme = replace(readme, "\r\n" => "\n")
    return replace(readme, r"<!-- BEGIN (\S+) -->.*?<!-- END \1 -->"s => s -> readme_block(match(r"<!-- BEGIN (\S+) -->", s)[1]))
end

if abspath(PROGRAM_FILE) == @__FILE__
    write(README_PATH, readme_with_examples(read(README_PATH, String)))
end
