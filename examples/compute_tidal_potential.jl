using Oceananigans
using Tidejinks
using Dates
import SPICE
using GLMakie

grid = LatitudeLongitudeGrid(size=(360, 180, 1), longitude=(0, 360), latitude=(-90, 90), z=(0, 1))
Φ = Field{Center, Center, Nothing}(grid)
t = DateTime(1993, 1, 1, 1)

kernel_meta_file = "kernels.txt"
Tidejinks.wrangle_spice_kernels(kernel_meta_file)
SPICE.furnsh(kernel_meta_file)

Tidejinks.compute_tidal_potential!(Φ, t)

heatmap(Φ)

