# # Visualize fields on a latitude-longitude grid.

# ## Install dependencies
#
# First let's make sure we have all required packages installed.

# ```julia
# using Pkg
# pkg"add Oceananigans, CairoMakie, Imaginocean, GeoMakie"
# ```

# Let's plot a field that lives on a latitude-longitude grid.

using Oceananigans

# First create a latitude-longitude grid

Nx, Ny, Nz = 80, 60, 2

grid = LatitudeLongitudeGrid(size = (Nx, Ny, Nz),
                             latitude = (-60, 60),
                             longitude = (-150, 20),
                             z = (-1, 0),
                             topology = (Bounded, Bounded, Bounded))

# Let's create a field. We choose here a field that lives on the ``y``-faces of the cells
# but any field would do.
# 
# We set the field value to ``\cos(3λ)^2 \sin(3φ)`` and see how that looks.

field = YFaceField(grid)

set!(field, (λ, φ, z) -> cosd(3λ)^2 * sind(3φ))

# We can visualize this field in 2D using a heatmap. Imaginocean.jl exports
# a method for `heatmap!` that works with Oceananigans.jl fields.

using CairoMakie, Imaginocean

kwargs = (colorrange = (-1, 1), colormap = :balance)

fig = Figure()
ax = Axis(fig[1, 1],
          xlabel="longitude [ᵒ]",
          ylabel="latitude [ᵒ]")

heatmap!(ax, field; kwargs...)

current_figure() # hide
fig


# We can do the same but with a `GeoAxis` provided by the GeoMakie.jl package
# that allows us to add coastlines easy or use various projections.

using GeoMakie

fig = Figure()
ax = GeoAxis(fig[1, 1],
             coastlines = true,
             lonlims = automatic)

heatmap!(ax, field; kwargs...)

current_figure() # hide
fig

# To make a 3D visualization on the sphere we first create a 3D axis and then
# use `heatsphere!` method from Imaginocean.jl.

fig = Figure()
ax = Axis3(fig[1, 1], aspect=(1, 1, 1), limits=((-1, 1), (-1, 1), (-1, 1)))

heatsphere!(ax, field; kwargs...)
hidedecorations!(ax) # hides the axes labels

current_figure() # hide
fig
