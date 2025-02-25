using ClimaOcean
using ClimaOcean.DataWrangling.JRA55: JRA55_field_time_series
using Oceananigans
using Tidejinks
using Dates
using GLMakie

backend = JRA55NetCDFBackend(41)
pa = JRA55_field_time_series(:sea_level_pressure)

times = pa.times
grid = pa.grid
Φ = Field{Center, Center, Nothing}(grid)

t = DateTime(1993, 1, 1, 1) + Second(times[1])

import SPICE
kernel_meta_file = "kernels.txt"
Tidejinks.wrangle_spice_kernels(kernel_meta_file)
SPICE.furnsh(kernel_meta_file)

Tidejinks.compute_tidal_potential!(Φ, t)

heatmap(Φ)
