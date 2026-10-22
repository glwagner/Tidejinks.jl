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

    for reference_frame in ("IAU_EARTH", "ITRF93"), date in unique(row[2:6] for row in rows)
        t = DateTime(parse.(Int, date)...)
        reference = [parse(Float64, row[9]) for row in rows if row[1] == reference_frame && row[2:6] == date]
        Tidejinks.compute_tidal_potential!(Φ, t, reference_frame)
        @test vec(interior(Φ)) ≈ reference atol=1e-9
    end
end


@testset "TidalPotential against the direct formula" begin
    start_date = DateTime(2015, 1, 1)
    dates = (start_date, start_date + Day(30))
    G_sun, G_moon = Tidejinks.gravitational_parameters()

    # Model times up to 2.6e7 s, which Float32 resolves to 2 s, so these times are exact in both precisions
    for FT in (Float64, Float32), epoch in (start_date, start_date - Day(270))
        grid = LatitudeLongitudeGrid(FT, size=(4, 3, 1), longitude=(0, 360), latitude=(-60, 60), z=(0, 1))
        Φ = Tidejinks.TidalPotential(grid, dates; epoch, reference_frame="IAU_EARTH")
        error = 0.0
        for date in start_date .+ Second.(2 .* rand(0:15 * 86400, 20)), λ in 0:30:330, φ in -90:30:90
            t = convert(FT, Dates.value(date - epoch) / 1000)
            X_sun, X_moon = Tidejinks.celestial_positions(date, "IAU_EARTH")
            exact = Tidejinks.compute_tidal_potential(deg2rad(λ), deg2rad(φ), X_sun, X_moon, G_sun, G_moon)
            error = max(error, abs(Φ(convert(FT, λ), convert(FT, φ), t) - exact))
        end
        @test error < (FT == Float64 ? 1e-7 : 1e-5)
    end
end
