using ClimaOcean
using ClimaOcean.DataWrangling.JRA55: JRA55_field_time_series
using Oceananigans
using Tidejinks
using Dates
using GLMakie
using SPICE

kernel_meta_file = "kernels.txt"
Tidejinks.wrangle_spice_kernels(kernel_meta_file)
SPICE.furnsh(kernel_meta_file)

backend = JRA55NetCDFBackend(41)
pa = JRA55_field_time_series(:sea_level_pressure; backend)
times = pa.times
grid = pa.grid

filename = "tidal_potential_jra55.jld2"
Φt = FieldTimeSeries{Center, Center, Nothing}(grid, times, backend=OnDisk(),
                                              path = filename, name="Φ")

Φ = Field{Center, Center, Nothing}(grid)
Nt = length(times)
for n = 1:Nt
    t = DateTime(1993, 1, 1, 1) + Second(times[n])
    @info "Computing tidal potential for $t..."
    Tidejinks.compute_tidal_potential!(Φ, t)
    set!(Φt, Φ, n)
end

