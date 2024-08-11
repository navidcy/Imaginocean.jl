var documenterSearchIndex = {"docs":
[{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"EditURL = \"../../../examples/ConformalCubedSphere_grid.jl\"","category":"page"},{"location":"literated/ConformalCubedSphere_grid/#Visualize-fields-on-a-latitude-longitude-grid.","page":"Conformal cubed sphere grid","title":"Visualize fields on a latitude-longitude grid.","text":"","category":"section"},{"location":"literated/ConformalCubedSphere_grid/#Install-dependencies","page":"Conformal cubed sphere grid","title":"Install dependencies","text":"","category":"section"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"First let's make sure we have all required packages installed.","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"using Pkg\npkg\"add Oceananigans, GLMakie, Imaginocean, GeoMakie\"","category":"page"},{"location":"literated/ConformalCubedSphere_grid/#Construct-a-test-bed-field","page":"Conformal cubed sphere grid","title":"Construct a test-bed field","text":"","category":"section"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"Let's plot a field that lives on a cubed sphere grid.","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"using Oceananigans","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"First create a conformal cubed sphere grid.","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"Nx = 30\nNy = 30\nNz = 1\n\nradius = 1\n\ngrid = ConformalCubedSphereGrid(; panel_size = (Nx, Ny, Nz),\n                                  z = (-1, 0),\n                                  radius)","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"Let's create a field. We choose a field that lives on the center of the cells.","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"We set the field values to something and see how that looks.","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"field = CenterField(grid)\n\nset!(field, (λ, φ, z) -> (sind(3λ) + 1/3 * sind(5λ)) * cosd(3φ)^2)","category":"page"},{"location":"literated/ConformalCubedSphere_grid/#2D-visualization","page":"Conformal cubed sphere grid","title":"2D visualization","text":"","category":"section"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"We can visualize this field in 2D using a heatmap. Imaginocean.jl has a method called heatlatlon! that plots a field that lives on a grid whose native coordinates are latitude and longitude.","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"using GLMakie, Imaginocean\n\nkwargs = (colorrange = (-1, 1), colormap = :balance)\n\nfig = Figure()\nax = Axis(fig[1, 1],\n          xlabel = \"longitude [ᵒ]\",\n          ylabel = \"latitude [ᵒ]\",\n          limits = ((-180, 180), (-90, 90)))\n\nheatlatlon!(ax, field, 1; kwargs...)\n\ncurrent_figure()","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"We can do the same but with a GeoAxis provided by the GeoMakie.jl package that allows us to easily add coastlines or also use various projections.","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"using GeoMakie\n\nfig = Figure()\nax = GeoAxis(fig[1, 1],\n             coastlines = true,\n             lonlims = automatic)\n\nheatlatlon!(ax, field, 1; kwargs...)\n\ncurrent_figure()","category":"page"},{"location":"literated/ConformalCubedSphere_grid/#3D-visualization-on-the-sphere","page":"Conformal cubed sphere grid","title":"3D visualization on the sphere","text":"","category":"section"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"To make a 3D visualization on the sphere we first create a 3D axis and then use heatsphere! method from Imaginocean.jl.","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"fig = Figure()\nax = Axis3(fig[1, 1],\n           aspect = (1, 1, 1),\n           limits = ((-1, 1), (-1, 1), (-1, 1)))\n\nheatsphere!(ax, field; kwargs...)\nhidedecorations!(ax) # hides the axes labels\n\ncurrent_figure()","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"","category":"page"},{"location":"literated/ConformalCubedSphere_grid/","page":"Conformal cubed sphere grid","title":"Conformal cubed sphere grid","text":"This page was generated using Literate.jl.","category":"page"},{"location":"library/function_index/#main-index","page":"Function index","title":"Index","text":"","category":"section"},{"location":"library/function_index/","page":"Function index","title":"Function index","text":"Pages = [\"public.md\", \"internals.md\", \"function_index.md\"]","category":"page"},{"location":"library/outline/#Library-Outline","page":"Contents","title":"Library Outline","text":"","category":"section"},{"location":"library/outline/","page":"Contents","title":"Contents","text":"Pages = [\"public.md\", \"internals.md\", \"function_index.md\"]","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"EditURL = \"../../../examples/LatitudeLongitude_grid.jl\"","category":"page"},{"location":"literated/LatitudeLongitude_grid/#Visualize-fields-on-a-latitude-longitude-grid.","page":"Latitude-Longitude grid","title":"Visualize fields on a latitude-longitude grid.","text":"","category":"section"},{"location":"literated/LatitudeLongitude_grid/#Install-dependencies","page":"Latitude-Longitude grid","title":"Install dependencies","text":"","category":"section"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"First let's make sure we have all required packages installed.","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"using Pkg\npkg\"add Oceananigans, GLMakie, Imaginocean, GeoMakie\"","category":"page"},{"location":"literated/LatitudeLongitude_grid/#Construct-a-test-bed-field","page":"Latitude-Longitude grid","title":"Construct a test-bed field","text":"","category":"section"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"Let's plot a field that lives on a latitude-longitude grid.","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"using Oceananigans","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"First create a latitude-longitude grid.","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"Nx, Ny, Nz = 180, 120, 2\n\ngrid = LatitudeLongitudeGrid(size = (Nx, Ny, Nz),\n                             latitude = (-60, 60),\n                             longitude = (-155, 25),\n                             z = (-1, 0),\n                             topology = (Bounded, Bounded, Bounded))","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"Let's create a field. We choose a field that lives on the faces of the cells but any field should do.","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"We set the field value to sin^2(3λ) sin(3φ) and see how that looks.","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"field = Field{Face, Face, Center}(grid)\n\nset!(field, (λ, φ, z) -> sind(3λ)^2 * sind(3φ))","category":"page"},{"location":"literated/LatitudeLongitude_grid/#2D-visualization","page":"Latitude-Longitude grid","title":"2D visualization","text":"","category":"section"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"We can visualize this field in 2D using a heatmap. Imaginocean.jl adds a method to heatmap! so that it works with Oceananigans' fields.","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"using GLMakie, Imaginocean\n\nkwargs = (colorrange = (-1, 1), colormap = :balance)\n\nfig = Figure()\nax = Axis(fig[1, 1],\n          xlabel = \"longitude [ᵒ]\",\n          ylabel = \"latitude [ᵒ]\",\n          limits = ((-180, 180), (-90, 90)))\n\nheatmap!(ax, field, 1; kwargs...)\n\ncurrent_figure()","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"We can do the same but with a GeoAxis provided by the GeoMakie.jl package that allows us to easily add coastlines or also use various projections.","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"using GeoMakie\n\nfig = Figure()\nax = GeoAxis(fig[1, 1],\n             coastlines = true,\n             lonlims = automatic)\n\nheatlatlon!(ax, field, 1; kwargs...)\n\ncurrent_figure()","category":"page"},{"location":"literated/LatitudeLongitude_grid/#3D-visualization-on-the-sphere","page":"Latitude-Longitude grid","title":"3D visualization on the sphere","text":"","category":"section"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"To make a 3D visualization on the sphere we first create a 3D axis and then use heatsphere! method from Imaginocean.jl.","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"fig = Figure()\nax = Axis3(fig[1, 1],\n           aspect = (1, 1, 1),\n           limits = ((-1, 1), (-1, 1), (-1, 1)))\n\nheatsphere!(ax, field; kwargs...)\nhidedecorations!(ax) # hides the axes labels\n\ncurrent_figure()","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"","category":"page"},{"location":"literated/LatitudeLongitude_grid/","page":"Latitude-Longitude grid","title":"Latitude-Longitude grid","text":"This page was generated using Literate.jl.","category":"page"},{"location":"library/internals/#Private-types-and-functions","page":"Private","title":"Private types and functions","text":"","category":"section"},{"location":"library/internals/","page":"Private","title":"Private","text":"Documentation for Imaginocean.jl's internal interface.","category":"page"},{"location":"library/internals/#Imaginocean","page":"Private","title":"Imaginocean","text":"","category":"section"},{"location":"library/internals/","page":"Private","title":"Private","text":"Modules = [Imaginocean]\nPublic = false","category":"page"},{"location":"library/internals/#Imaginocean.get_cartesian_nodes_and_vertices-Tuple{Union{Oceananigans.Grids.LatitudeLongitudeGrid, Oceananigans.Grids.OrthogonalSphericalShellGrid}, Any, Any, Any}","page":"Private","title":"Imaginocean.get_cartesian_nodes_and_vertices","text":"get_cartesian_nodes_and_vertices(grid::Union{LatitudeLongitudeGrid, OrthogonalSphericalShellGrid}, ℓx, ℓy, ℓz)\n\nReturn the cartesian coordinates of the horizontal nodes of the grid at locations ℓx, ℓy, and ℓz on the unit sphere and also the corresponding coordinates of the four vertices that determine the cell surrounding each node.\n\nSee get_lat_lon_nodes_and_vertices.\n\n\n\n\n\n","category":"method"},{"location":"library/internals/#Imaginocean.get_lat_lon_nodes_and_vertices-NTuple{4, Any}","page":"Private","title":"Imaginocean.get_lat_lon_nodes_and_vertices","text":"get_lat_lon_nodes_and_vertices(grid, ℓx, ℓy, ℓzs; lower_limit=-180)\n\nReturn the latitude-longitude coordinates of the horizontal nodes of the grid at locations ℓx, ℓy, and ℓz and also the coordinates of the four vertices that determine the cell surrounding each node.\n\nSee get_longitude_vertices and get_latitude_vertices.\n\n\n\n\n\n","category":"method"},{"location":"library/internals/#Imaginocean.get_latitude_vertices-Tuple{Any, Any, Any, Union{Oceananigans.Grids.LatitudeLongitudeGrid, Oceananigans.Grids.OrthogonalSphericalShellGrid}, Any, Any, Any}","page":"Private","title":"Imaginocean.get_latitude_vertices","text":"get_latitude_vertices(i, j, k, grid::Union{LatitudeLongitudeGrid, OrthogonalSphericalShellGrid}, ℓx, ℓy, ℓz)\n\nReturn the latitudes that correspond to the four vertices of cell i, j, k at location (ℓx, ℓy, ℓz). The first vertex is the cell's Southern-Western oneλand the rest follow in counter-clockwise order.\n\n\n\n\n\n","category":"method"},{"location":"library/internals/#Imaginocean.get_longitude_vertices-Tuple{Any, Any, Any, Union{Oceananigans.Grids.LatitudeLongitudeGrid, Oceananigans.Grids.OrthogonalSphericalShellGrid}, Any, Any, Any}","page":"Private","title":"Imaginocean.get_longitude_vertices","text":"get_longitude_vertices(i, j, k, grid::Union{LatitudeLongitudeGrid, OrthogonalSphericalShellGrid}, ℓx, ℓy, ℓz)\n\nReturn the longitudes that correspond to the four vertices of cell i, j, k at location (ℓx, ℓy, ℓz). The first vertex is the cell's Southern-Western one and the rest follow in counter-clockwise order.\n\n\n\n\n\n","category":"method"},{"location":"library/internals/#Imaginocean.lat_lon_to_cartesian-Tuple{Any, Any}","page":"Private","title":"Imaginocean.lat_lon_to_cartesian","text":"lat_lon_to_cartesian(longitude, latitude; radius=1)\n\nConvert (λ φ) = (longitudelatitude) coordinates (in degrees) to cartesian coordinates (x y z) on a sphere with radius, R, i.e.,\n\nbeginaligned\nx = R cos(λ) cos(φ) \ny = R sin(λ) cos(φ) \nz = R sin(φ)\nendaligned\n\n\n\n\n\n","category":"method"},{"location":"library/internals/#Imaginocean.lat_lon_to_x-Tuple{Any, Any}","page":"Private","title":"Imaginocean.lat_lon_to_x","text":"lat_lon_to_x(longitude, latitude)\n\nConvert (longitude, latitude) coordinates (in degrees) to cartesian x on the unit sphere.\n\n\n\n\n\n","category":"method"},{"location":"library/internals/#Imaginocean.lat_lon_to_y-Tuple{Any, Any}","page":"Private","title":"Imaginocean.lat_lon_to_y","text":"lat_lon_to_y(longitude, latitude)\n\nConvert (longitude, latitude) coordinates (in degrees) to cartesian y on the unit sphere.\n\n\n\n\n\n","category":"method"},{"location":"library/internals/#Imaginocean.lat_lon_to_z-Tuple{Any, Any}","page":"Private","title":"Imaginocean.lat_lon_to_z","text":"lat_lon_to_z(longitude, latitude)\n\nConvert (longitude, latitude) coordinates (in degrees) to cartesian z on the unit sphere.\n\n\n\n\n\n","category":"method"},{"location":"library/internals/#Imaginocean.longitude_domain-Tuple{Any}","page":"Private","title":"Imaginocean.longitude_domain","text":"longitude_domain(longitude; lower_limit = -180)\n\nBring longitude to domain [lower_limit, lower_limit + 360] (in degrees). By default, lower_limit = -180 implying longitude domain -180 180.\n\nExamples\n\njulia> using Imaginocean: longitude_domain\n\njulia> longitude_domain(400)\n40\n\njulia> longitude_domain(-50)\n-50\n\njulia> longitude_domain(-50; lower_limit=0)\n310\n\n\n\n\n\n","category":"method"},{"location":"library/internals/#MakieCore.convert_arguments-Tuple{Type{<:AbstractPlot}, Oceananigans.Fields.Field, Int64}","page":"Private","title":"MakieCore.convert_arguments","text":"Makie.convert_arguments(P::Type{<:AbstractPlot}, field::Field, k_index::Int)\n\nConvert an Oceananigans.jl Field with non-flat horizontal dimensions at vertical index k_index to arguments that can be plotted as a surface in Makie.jl.\n\n\n\n\n\n","category":"method"},{"location":"#Imaginocean.jl","page":"Home","title":"Imaginocean.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Visualization package for Oceananigans.jl fields using the plotting package Makie.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Makie.jl comes with a few backends. In the documented examples we use GLMakie backend.","category":"page"},{"location":"library/public/#Public-Documentation","page":"Public","title":"Public Documentation","text":"","category":"section"},{"location":"library/public/","page":"Public","title":"Public","text":"Documentation for Imaginocean.jl's public interface.","category":"page"},{"location":"library/public/","page":"Public","title":"Public","text":"See the Internals section of the manual for internal package docs covering all submodules.","category":"page"},{"location":"library/public/#Imaginocean","page":"Public","title":"Imaginocean","text":"","category":"section"},{"location":"library/public/","page":"Public","title":"Public","text":"Modules = [Imaginocean]\nPrivate = false","category":"page"},{"location":"library/public/#Imaginocean.heatsphere!","page":"Public","title":"Imaginocean.heatsphere!","text":"heatsphere!(axis::Axis3, field::Field, k_index=1; kwargs...)\n\nA heatmap of an Oceananigans.jl Field on the sphere at vertical index k_index.\n\nArguments\n\naxis :: Makie.Axis3: a 3D axis.\nfield :: Oceananigans.Field: an Oceananigans.jl field with non-flat horizontal dimensions.\nk_index :: Int: The integer corresponding to the vertical index of the field to visualize; default: 1.\n\nAccepts all keyword arguments for Makie.mesh! method.\n\n\n\n\n\n","category":"function"}]
}
