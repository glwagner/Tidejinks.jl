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

@testset "Computing tidal potentials..." begin
    grid = LatitudeLongitudeGrid(size=(90, 45, 1), longitude=(0, 360), latitude=(-90, 90), z=(0, 1))
    Φ = Field{Center, Center, Nothing}(grid)
    t = DateTime(1993, 1, 1, 1)
    Tidejinks.compute_tidal_potential!(Φ, t)
    @test maximum(Φ) ≈ 1.6422346052783796
    @test minimum(Φ) ≈ -1.5737500023619164
end


