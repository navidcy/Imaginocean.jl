# Imaginocean.jl

Visualization package for Oceananigans.jl fields using [Makie.jl](https://docs.makie.org/stable/) for plotting.

Makie comes with a few [backends](https://docs.makie.org/stable/#makie_ecosystem). In the examples
we use [CairoMakie](https://docs.makie.org/stable/documentation/backends/cairomakie/) since this backend
works well on headless devices, that is, devices without monitor. The documentation is automatically
built via GitHub actions and thus CairoMakie backend is necessary.

Users that using devices with a monitor might want to change to
[GLMakie](https://docs.makie.org/stable/documentation/backends/glmakie/)
that displays figures in an interactive window.
