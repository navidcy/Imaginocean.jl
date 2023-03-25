# Imaginocean.jl

Visualization package for [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/)
fields using the plotting package [Makie.jl](https://docs.makie.org/stable/).

Makie.jl comes with a few [backends](https://docs.makie.org/stable/#makie_ecosystem).
In the documented examples we use [CairoMakie](https://docs.makie.org/stable/documentation/backends/cairomakie/)
since this backend works well on headless devices, that is, devices without monitor. Since the
package's documentation is built automatically via GitHub actions using CairoMakie backend is
necessary.

Users on devices with a monitor might want to change to
[GLMakie](https://docs.makie.org/stable/documentation/backends/glmakie/) backend
which will show the figures in an interactive window.
