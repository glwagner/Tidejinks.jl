using Oceananigans
using Tidejinks
using Dates
using GLMakie
using SPICE

kernel_meta_file = "kernels.txt"
Tidejinks.wrangle_spice_kernels(kernel_meta_file)
SPICE.furnsh(kernel_meta_file)

grid = LatitudeLongitudeGrid(size=(360, 180, 1), longitude=(0, 360), latitude=(-90, 90), z=(0, 1))
Φ = Field{Center, Center, Nothing}(grid)

t₀ = DateTime(1993, 1, 1, 1)
t₁ = t₀ + Day(60)
dt = Hour(3)
times = t₀:dt:t₁
Nt = length(times)

Φt = []
for t in times
    @info "Computing tidal potential at t = $t..."
    Tidejinks.compute_tidal_potential!(Φ, t)
    push!(Φt, deepcopy(Φ))
end

fig = Figure()
ax = Axis(fig[2, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)", aspect=2)
n = Observable(1)
Φn = @lift Φt[$n]
hm = heatmap!(ax, Φn)
Colorbar(fig[3, 1], hm, label="Tidal potential at Earth's surface (m² s⁻²)", vertical=false)
titlestr = @lift string("Time: ", times[$n])
Label(fig[1, 1], titlestr, tellwidth=false)

record(fig, "sixty_day_tidal_potential.mp4", 1:Nt, framerate=12) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end

