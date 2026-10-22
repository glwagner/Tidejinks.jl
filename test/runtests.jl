using Test
using Tidejinks
using Dates
using Oceananigans
import SPICE

@testset "Test kernel wrangling..." begin
    kernel_meta_file = "kernels.txt"
    Tidejinks.wrangle_spice_kernels(kernel_meta_file)
    SPICE.furnsh(kernel_meta_file)
    @test isfile(kernel_meta_file)
end

@testset "Tidal potential against MITgcm's tides_calc.F" begin
    rows = [split(line) for line in eachline(joinpath(@__DIR__, "tides_calc_reference.dat")) if !startswith(line, "#")]
    grid = LatitudeLongitudeGrid(size=(4, 3, 1), longitude=(0, 360), latitude=(-60, 60), z=(0, 1))
    Φ = Field{Center, Center, Nothing}(grid)

    for frame in ("IAU_EARTH", "ITRF93"), date in unique(row[2:6] for row in rows)
        t = DateTime(parse.(Int, date)...)
        reference = [parse(Float64, row[9]) for row in rows if row[1] == frame && row[2:6] == date]
        Tidejinks.compute_tidal_potential!(Φ, t, frame)
        @test vec(interior(Φ)) ≈ reference atol=1e-9
    end
end

