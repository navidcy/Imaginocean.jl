using Oceananigans

using Visualizocean, GLMakie, GeoMakie

GLMakie.activate!()

Nx, Ny, Nz = 80, 60, 2

grid = LatitudeLongitudeGrid(size = (Nx, Ny, Nz),
                             latitude = (-60, 60),
                             longitude = (-150, 20),
                             z = (-1, 0),
                             topology = (Bounded, Bounded, Bounded),
                             halo = (2, 2, 2))

field = YFaceField(grid)

set!(field, (λ, φ, z) -> cosd(3λ)^2 * sind(3φ))
colorrange = (-1, 1)
colormap = :balance

fig = Figure()
ax = Axis3(fig[1, 1], aspect=(1, 1, 1), limits=((-1, 1), (-1, 1), (-1, 1)))

heatsphere!(ax, field; colorrange, colormap)

save("lat_lon_on_sphere.png", fig)
fig

fig = Figure()
ax = Axis(fig[1, 1], xlabel="longitude [ᵒ]", ylabel="latitude [ᵒ]")

heatlatlon!(ax, field; colorrange, colormap)

save("lat_lon_2D.png", fig)
fig


fig = Figure()
ax = GeoAxis(fig[1, 1],
             coastlines = true,
             lonlims = automatic)

heatlatlon!(ax, field; colorrange, colormap)

save("lat_lon_GeoMakie.png", fig)
fig
