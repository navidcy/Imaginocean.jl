# # Visualize fields on a latitude-longitude grid.

# ### Install dependencies
#
# First let's make sure we have all required packages installed.

# ```julia
# using Pkg
# pkg"add Oceananigans, CairoMakie, Imaginocean, GeoMakie"
# ```

# ### Construct a test-bed field

# Let's plot a field that lives on a cubed sphere grid.

using Oceananigans

# First create a conformal cubed sphere grid.

Nx = 30
Ny = 30
Nz = 1

radius = 1

grid = ConformalCubedSphereGrid(; panel_size = (Nx, Ny, Nz),
                                  z = (-1, 0),
                                  radius)

# Let's create a field. We choose a field that lives on the center of the cells.
#
# We set the field values to something and see how that looks.

field = CenterField(grid)

set!(field, (λ, φ, z) -> (sind(3λ) + 1/3 * sind(5λ)) * cosd(3φ)^2)

# ### 2D visualization

# We can visualize this field in 2D using a heatmap. Imaginocean.jl has a method
# called `heatlatlon!` which plots a field that lives on a grid whose native
# coordinates are latitude-longitude.

using CairoMakie, Imaginocean

kwargs = (colorrange = (-1, 1), colormap = :balance)

fig = Figure()
ax = Axis(fig[1, 1],
          xlabel = "longitude [ᵒ]",
          ylabel = "latitude [ᵒ]",
          limits = ((-180, 180), (-90, 90)))

heatlatlon!(ax, field, 1; kwargs...)

current_figure() #hide
fig

# We can do the same but with a `GeoAxis` provided by the GeoMakie.jl package
# that allows us to easily add coastlines or also use various projections.

using GeoMakie

fig = Figure()
ax = GeoAxis(fig[1, 1],
             coastlines = true,
             lonlims = automatic)

heatlatlon!(ax, field, 1; kwargs...)

current_figure() # hide
fig

# ### 3D visualization on the sphere

# To make a 3D visualization on the sphere we first create a 3D axis and then
# use `heatsphere!` method from Imaginocean.jl.

fig = Figure()
ax = Axis3(fig[1, 1],
           aspect = (1, 1, 1),
           limits = ((-1, 1), (-1, 1), (-1, 1)))

heatsphere!(ax, field; kwargs...)
hidedecorations!(ax) # hides the axes labels

current_figure() #hide
fig
