module Imaginocean

export heatsphere!

using Makie
using Oceananigans
using Oceananigans.Grids: xnode, ynode, total_length, topology

"""
    lat_lon_to_cartesian(longitude, latitude; radius=1)

Convert ``(λ, φ) =(```longitude```, ```latitude```)`` coordinates (in degrees) to
cartesian coordinates ``(x, y, z)`` on a sphere with `radius`, ``R``, i.e.,

```math
\\begin{aligned}
x &= R \\cos(λ) \\cos(φ), \\\\
y &= R \\sin(λ) \\cos(φ), \\\\
z &= R \\sin(φ).
\\end{aligned}
```
"""
lat_lon_to_cartesian(longitude, latitude; radius=1) = (radius * lat_lon_to_x(longitude, latitude),
                                                       radius * lat_lon_to_y(longitude, latitude),
                                                       radius * lat_lon_to_z(longitude, latitude))

"""
    lat_lon_to_x(longitude, latitude)

Convert `(longitude, latitude)` coordinates (in degrees) to cartesian `x` on the unit sphere.
"""
lat_lon_to_x(longitude, latitude) = cosd(longitude) * cosd(latitude)

"""
    lat_lon_to_y(longitude, latitude)

Convert `(longitude, latitude)` coordinates (in degrees) to cartesian `y` on the unit sphere.
"""
lat_lon_to_y(longitude, latitude) = sind(longitude) * cosd(latitude)

"""
    lat_lon_to_z(longitude, latitude)

Convert `(longitude, latitude)` coordinates (in degrees) to cartesian `z` on the unit sphere.
"""
lat_lon_to_z(longitude, latitude) = sind(latitude)

"""
    longitude_domain(longitude; lower_limit = -180)

Bring `longitude` to domain `[lower_limit, lower_limit+360]` (in degrees).
By default, `lower_limit = -180` implying longitude domain ``[-180, 180]``.

Examples
========

```jldoctest
julia> using Imaginocean: longitude_domain

julia> longitude_domain(400)
40

julia> longitude_domain(-50)
-50

julia> longitude_domain(-50; lower_limit=0)
310
```
"""
longitude_domain(longitude; lower_limit = -180) = mod.(longitude .+ lower_limit + 360, 360) .+ lower_limit

flip_location(::Center) = Face()
flip_location(::Face) = Center()

"""
    get_longitude_vertices(i, j, k, grid::Union{LatitudeLongitudeGrid, OrthogonalSphericalShellGrid}, ℓx, ℓy, ℓz)

Return the longitudes that correspond to the four vertices of cell `i, j, k` at
location `(ℓx, ℓy, ℓz)`. The first vertex is the cell's Southern-Western one
and the rest follow in counter-clockwise order.
"""
function get_longitude_vertices(i, j, k, grid::Union{LatitudeLongitudeGrid, OrthogonalSphericalShellGrid}, ℓx, ℓy, ℓz)

    if ℓx == Center()
        i₀ = i
    elseif ℓx == Face()
        i₀ = i-1
    end

    if ℓy == Center()
        j₀ = j
    elseif ℓy == Face()
        j₀ = j-1
    end

    λ_vertex₁ = xnode( i₀,   j₀,  k, grid, flip_location(ℓx), flip_location(ℓy), ℓz)
    λ_vertex₂ = xnode(i₀+1,  j₀,  k, grid, flip_location(ℓx), flip_location(ℓy), ℓz)
    λ_vertex₃ = xnode(i₀+1, j₀+1, k, grid, flip_location(ℓx), flip_location(ℓy), ℓz)
    λ_vertex₄ = xnode( i₀,  j₀+1, k, grid, flip_location(ℓx), flip_location(ℓy), ℓz)

    return [λ_vertex₁; λ_vertex₂; λ_vertex₃; λ_vertex₄]
end

"""
    get_latitude_vertices(i, j, k, grid::Union{LatitudeLongitudeGrid, OrthogonalSphericalShellGrid}, ℓx, ℓy, ℓz)

Return the latitudes that correspond to the four vertices of cell `i, j, k` at
location `(ℓx, ℓy, ℓz)`. The first vertex is the cell's Southern-Western one
and the rest follow in counter-clockwise order.
"""
function get_latitude_vertices(i, j, k, grid::Union{LatitudeLongitudeGrid, OrthogonalSphericalShellGrid}, ℓx, ℓy, ℓz)

    if ℓx == Center()
        i₀ = i
    elseif ℓx == Face()
        i₀ = i-1
    end

    if ℓy == Center()
        j₀ = j
    elseif ℓy == Face()
        j₀ = j-1
    end

    φ_vertex₁ = ynode( i₀,   j₀,  k, grid, flip_location(ℓx), flip_location(ℓy), ℓz)
    φ_vertex₂ = ynode(i₀+1,  j₀,  k, grid, flip_location(ℓx), flip_location(ℓy), ℓz)
    φ_vertex₃ = ynode(i₀+1, j₀+1, k, grid, flip_location(ℓx), flip_location(ℓy), ℓz)
    φ_vertex₄ = ynode( i₀,  j₀+1, k, grid, flip_location(ℓx), flip_location(ℓy), ℓz)

    return [φ_vertex₁; φ_vertex₂; φ_vertex₃; φ_vertex₄]
end

longitude_in_same_window(λ₁, λ₂) = mod(λ₁ - λ₂ + 180, 360) + λ₂ - 180

"""
    get_lat_lon_nodes_and_vertices(grid, ℓx, ℓy, ℓzs; lower_limit=-180)

Return the latitude-longitude coordinates of the horizontal nodes of the
`grid` at locations `ℓx`, `ℓy`, and `ℓz` and also the coordinates of the four
vertices that determine the cell surrounding each node.

See [`get_longitude_vertices`](@ref) and [`get_latitude_vertices`](@ref).
"""
function get_lat_lon_nodes_and_vertices(grid, ℓx, ℓy, ℓz)

    TX, TY, _ = topology(grid)

    λ = zeros(eltype(grid), total_length(ℓx, TX(), grid.Nx, 0), total_length(ℓy, TY(), grid.Ny, 0))
    φ = zeros(eltype(grid), total_length(ℓx, TX(), grid.Nx, 0), total_length(ℓy, TY(), grid.Ny, 0))

    for j in axes(λ, 2), i in axes(λ, 1)
        λ[i, j] = xnode(i, j, 1, grid, ℓx, ℓy, ℓz)
        φ[i, j] = ynode(i, j, 1, grid, ℓx, ℓy, ℓz)
    end

    λvertices = zeros(4, size(λ)...)
    φvertices = zeros(4, size(φ)...)

    for j in axes(λ, 2), i in axes(λ, 1)
        λvertices[:, i, j] = get_longitude_vertices(i, j, 1, grid, ℓx, ℓy, ℓz)
        φvertices[:, i, j] =  get_latitude_vertices(i, j, 1, grid, ℓx, ℓy, ℓz)
    end
    
    # ensure λ ∈ [-180, 180]
    @. λ = longitude_domain(λ)

    # ensure all vertices have longitudes in the same domain as λ
    λvertices = longitude_domain.(λvertices .- reshape(λ, (1, size(λ)...))) .+ reshape(λ, (1, size(λ)...))

    return (λ, φ), (λvertices, φvertices)
end

"""
    get_cartesian_nodes_and_vertices(grid::Union{LatitudeLongitudeGrid, OrthogonalSphericalShellGrid}, ℓx, ℓy, ℓz)

Return the cartesian coordinates of the horizontal nodes of the `grid`
at locations `ℓx`, `ℓy`, and `ℓz` on the unit sphere and also the corresponding
coordinates of the four vertices that determine the cell surrounding each node.

See [`get_lat_lon_nodes_and_vertices`](@ref).
"""
function get_cartesian_nodes_and_vertices(grid::Union{LatitudeLongitudeGrid, OrthogonalSphericalShellGrid}, ℓx, ℓy, ℓz)

    (λ, φ), (λvertices, φvertices) = get_lat_lon_nodes_and_vertices(grid, ℓx, ℓy, ℓz)

    x = similar(λ)
    y = similar(λ)
    z = similar(λ)

    xvertices = similar(λvertices)
    yvertices = similar(λvertices)
    zvertices = similar(λvertices)

    for j in axes(λ, 2), i in axes(λ, 1)
        x[i, j] = lat_lon_to_x(λ[i, j], φ[i, j])
        y[i, j] = lat_lon_to_y(λ[i, j], φ[i, j])
        z[i, j] = lat_lon_to_z(λ[i, j], φ[i, j])

        for vertex in 1:4
            xvertices[vertex, i, j] = lat_lon_to_x(λvertices[vertex, i, j], φvertices[vertex, i, j])
            yvertices[vertex, i, j] = lat_lon_to_y(λvertices[vertex, i, j], φvertices[vertex, i, j])
            zvertices[vertex, i, j] = lat_lon_to_z(λvertices[vertex, i, j], φvertices[vertex, i, j])
        end
    end

    return (x, y, z), (xvertices, yvertices, zvertices)
end

"""
    heatsphere!(axis::Axis3, field::Field, k_index=1; kwargs...)

A heatmap of an Oceananigans.jl `Field` on the sphere at vertical index `k_index`.

Arguments
=========
* `axis :: Makie.Axis3`: a 3D axis.
* `field :: Oceananigans.Field`: an Oceananigans.jl field with non-flat horizontal dimensions.
* `k_index :: Int`: The integer corresponding to the vertical index of the `field` to visualize; default: 1.

Accepts all keyword arguments for `Makie.mesh!` method.
"""
function heatsphere!(axis::Axis3, field::Field, k_index=1; kwargs...)
    LX, LY, LZ = location(field)
    grid = field.grid

    _, (xvertices, yvertices, zvertices) = get_cartesian_nodes_and_vertices(grid, LX(), LY(), LZ())

    quad_points3 = vcat([Point3.(xvertices[:, i, j], yvertices[:, i, j], zvertices[:, i, j]) for i in axes(xvertices, 2), j in axes(xvertices, 3)]...)
    quad_faces = vcat([begin; j = (i-1) * 4 + 1; [j j+1  j+2; j+2 j+3 j]; end for i in 1:length(quad_points3)÷4]...)

    field_2D = interior(field, :, :, k_index)

    colors_per_point = vcat(fill.(vec(field_2D), 4)...)

    mesh!(axis, quad_points3, quad_faces; color = colors_per_point, shading = false, interpolate=false, kwargs...)

    return axis
end

"""
    Makie.convert_arguments(P::SurfaceLike, field::Field, k_index::Int)

Convert an Oceananigans.jl `Field` with non-flat horizontal dimensions at vertical
index `k_index` to arguments that can be plotted as a `SurfaceLike` type in Makie.jl.
"""
function Makie.convert_arguments(P::SurfaceLike, field::Field, k_index::Int; kwargs...)
    LX, LY, LZ = location(field)
    grid = field.grid

    _, (λvertices, φvertices) = get_lat_lon_nodes_and_vertices(grid, LX(), LY(), LZ())

    quad_points = vcat([Point2f.(λvertices[:, i, j], φvertices[:, i, j]) for i in axes(λvertices, 2), j in axes(λvertices, 3)]...)

    longitudes = unique([quad_points[i][1] for i ∈ eachindex(quad_points)])
    latitudes  = unique([quad_points[i][2] for i ∈ eachindex(quad_points)])

    field_2D = interior(field, :, :, k_index)

    return convert_arguments(P, longitudes, latitudes, field_2D)
end

end # module Imaginocean
