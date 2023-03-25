pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..")) # add Imaginocean to environment stack

using
  Documenter,
  Literate,
  CairoMakie,  # so that Literate.jl does not capture precompilation output or warnings
  Glob,
  Imaginocean

#####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

to_be_literated = [
    "LatitudeLongitude_grid.jl"
]

for file in to_be_literated
    filepath = joinpath(EXAMPLES_DIR, file)
    Literate.markdown(filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor())
end

example_pages = [
    "Latitude-Longitude grid" => "literated/LatitudeLongitude_grid.md",
]

#####
##### Build and deploy docs
#####

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://clima.github.io/ClimaOceanDocumentation/dev/",
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

@info "Cleaning up temporary .jld2 and .nc files created by doctests..."

for file in vcat(glob("docs/*.jld2"), glob("docs/*.nc"))
    rm(file)
end

deploydocs(        repo = "github.com/navidcy/Imaginocean.jl.git",
                versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
              forcepush = true,
              devbranch = "main",
            push_preview = true)
