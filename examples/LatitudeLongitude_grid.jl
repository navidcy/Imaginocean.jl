# Let's plot a field that lives on a latitude-longitude grid.

using Oceananigans
using Visualizocean
using CairoMakie, GeoMakie

Nx, Ny, Nz = 80, 60, 2

grid = LatitudeLongitudeGrid(size = (Nx, Ny, Nz),
                             latitude = (-60, 60),
                             longitude = (-150, 20),
                             z = (-1, 0),
                             topology = (Bounded, Bounded, Bounded))

field = YFaceField(grid)

set!(field, (λ, φ, z) -> cosd(3λ)^2 * sind(3φ))

colorrange = (-1, 1)
colormap = :balance
nothing # hide

# To plot on the sphere we first create a figure with a 3D axis.

fig = Figure()
ax = Axis3(fig[1, 1], aspect=(1, 1, 1), limits=((-1, 1), (-1, 1), (-1, 1)))

heatsphere!(ax, field; colorrange, colormap)

current_figure() # hide
fig

# Now let's make a 2D plot

fig = Figure()
ax = Axis(fig[1, 1], xlabel="longitude [ᵒ]", ylabel="latitude [ᵒ]")

heatlatlon!(ax, field; colorrange, colormap)

current_figure() # hide
fig

# and another one

fig = Figure()
ax = GeoAxis(fig[1, 1],
             coastlines = true,
             lonlims = automatic)

heatlatlon!(ax, field; colorrange, colormap)

current_figure() # hide
fig
