pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..")) # add Imaginocean to environment stack

using
  Documenter,
  Literate,
  CairoMakie,  # so that Literate.jl does not capture precompilation output or warnings
  Imaginocean

CairoMakie.activate!(type = "svg")

#####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

to_be_literated = [
    "LatitudeLongitude_grid.jl",
    "ConformalCubedSphere_grid.jl"
]

for file in to_be_literated
    filepath = joinpath(EXAMPLES_DIR, file)
    Literate.markdown(filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor())
end

example_pages = [
    "Latitude-Longitude grid" => "literated/LatitudeLongitude_grid.md",
    "Conformal cubed sphere grid" => "literated/ConformalCubedSphere_grid.md",
]

#####
##### Build and deploy docs
#####

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://navidcy.github.io/Imaginocean.jl/dev/",
)

pages = [
    "Home" => "index.md",
    "Examples" => example_pages,
    "Library" => [ 
        "Contents"       => "library/outline.md",
        "Public"         => "library/public.md",
        "Private"        => "library/internals.md",
        "Function index" => "library/function_index.md",
        ],
]

makedocs(
   sitename = "Imaginocean.jl",
    modules = [Imaginocean],
     format = format,
      pages = pages,
    doctest = true,
     strict = false,
      clean = true,
  checkdocs = :exports
)

@info "Clean up temporary .jld2 and .nc output created by doctests or literated examples..."

"""
    recursive_find(directory, pattern)

Return list of filepaths within `directory` that contains the `pattern::Regex`.
"""
recursive_find(directory, pattern) =
    mapreduce(vcat, walkdir(directory)) do (root, dirs, files)
        joinpath.(root, filter(contains(pattern), files))
    end

files = []
for pattern in [r"\.jld2", r"\.nc"]
    global files = vcat(files, recursive_find(@__DIR__, pattern))
end

for file in files
    rm(file)
end

deploydocs(       repo = "github.com/navidcy/Imaginocean.jl.git",
               versions = ["dev" => "dev", "stable" => "v^", "v#.#.#"],
             forcepush = true,
             devbranch = "main",
           push_preview = true)
