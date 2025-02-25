module Tidejinks

using Scratch
using Downloads
using ClimaOcean.DataWrangling: download_progress

urls = Dict(
    "de440.bsp" => "https://www.dropbox.com/scl/fi/rncpxkcad8fmr2oxboq0k/" *
                   "de440.bsp?rlkey=5ao2velhbqvb28mdzjs29pol2&st=sgya8xqa&dl=0",

    "earth_000101_220503_220207.bpc" => "https://www.dropbox.com/scl/fi/f6pedw7fkulba152st2sc/" *
                                        "earth_000101_220503_220207.bpc?rlkey=dm4k6u0xchtdu82fh9b40mqn5&st=lbc569eg&dl=0",

    "earth_720101_070426.bpc" => "https://www.dropbox.com/scl/fi/oh4jd57xy46ir3nq4lsnu/" *
                                 "earth_720101_070426.bpc?rlkey=v23cjs3mn5y1j2n0n6xq8fqv0&st=lyfz74d7&dl=0",

    "earth_assoc_itrf93.tf" => "https://www.dropbox.com/scl/fi/3qb0hjc2piymn2ouj6xvm/" * 
                               "earth_assoc_itrf93.tf?rlkey=0lrzv9tfgyha7pmxwxlf2cyw5&st=mrq7cxwl&dl=0",

    "gm_de431.tpc" => "https://www.dropbox.com/scl/fi/ks439kq6u5n0qe3kvxgr4/" *
                      "gm_de431.tpc?rlkey=tovdsv0j27pguuy2hdyuf1gpa&st=fri24ium&dl=0",

    "naif0012.tls" => "https://www.dropbox.com/scl/fi/btiewfxq46wqinxxdmu8k/" *
                      "naif0012.tls?rlkey=eep0924pqycbbnsib88al8mem&st=j9qba8zm&dl=0",

    "pck00010.tpc" => "https://www.dropbox.com/scl/fi/tqgp7jg4avm6ulyrwa0zh/" *
                      "pck00010.tpc?rlkey=johcphes1h8a2nkxml9euth15&st=vcv77231&dl=0",

    "latest_leapseconds.tls" => "https://www.dropbox.com/scl/fi/btcttj688gxuxki90mor6/" * 
                                "latest_leapseconds.tls?rlkey=66xhz970q7tzrtvtrmxto20ie&st=hwmbtfks&dl=0",
)

#spice_cache::String = ""
#meta_kernel_file::String = ""

spice_cache = "SPICE" #@get_scratch!("SPICE")
!ispath(spice_cache) && mkdir(spice_cache)
meta_kernel_file = joinpath(spice_cache, "meta_kernel.txt")

for (path, url) in urls
    filepath = joinpath(spice_cache, path)
    @info "Downloading file to $filepath..."
    Downloads.download(url, filepath; progress=download_progress)
end

meta_kernel = """
\\begindata

KERNELS_TO_LOAD=(
    '$spice_cache/naif0012.tls',
    '$spice_cache/latest_leapseconds.tls',
    '$spice_cache/earth_assoc_itrf93.tf',
    '$spice_cache/de440.bsp',
    '$spice_cache/pck00010.tpc',
    '$spice_cache/gm_de431.tpc',
    '$spice_cache/earth_000101_220503_220207.bpc'
    '$spice_cache/earth_720101_070426.bpc'
)

\\begintext
"""

if !isfile(meta_kernel_file)
    open(meta_kernel_file, "w") do file
        write(file, meta_kernel)
    end
end

using SPICE
@show spice_cache
@show meta_kernel_file
furnsh(meta_kernel_file)

# Get gravitational parameters (in m)
const G_sun  = first(bodvrd("SUN", "GM", 1)) * 1e9
const G_moon = first(bodvrd("MOON", "GM", 1)) * 1e9

# Constants
const EARTH_RADIUS = 6371e3
const LOVE_H2 = 0.61
const LOVE_K2 = 0.30

using Oceananigans
using Dates
using LegendrePolynomials

"""
    get_body_position(body, time)

Get longitude, latitude, and distance of celestial body.
"""
function get_body_position(body, time)
    time_str = Dates.format(time, "yyyy-mm-dd HH:MM:SS") * " UTC"
    ephemeris_time = str2et(time_str)
    position_rectangular, _ = spkpos(body, ephemeris_time, "ITRF93", "NONE", "EARTH")
    position_geodetic = recgeo(position_rectangular, EARTH_RADIUS/1e3, 0)
    position_radial = recrad(position_rectangular)

    @inbounds begin
        λ = position_geodetic[1]
        φ = position_geodetic[2]
        r = position_radial[1] * 1e3
    end

    return λ, φ, r
end

"""
    calculate_zenith_cosine(λ₁, φ₁, λ₂, φ₂)

Calculate cosine of zenith angle using spherical law of cosines.
"""
calculate_zenith_cosine(λ₁, φ₁, λ₂, φ₂) =
    sind(φ₁)*sind(φ₂) + cosd(φ₁)*cosd(φ₂)*cosd(λ₂ - λ₁)

"""
    compute_potential(longitude, latitude, time; density=1)

Compute tidal potential at specified time for given grid coordinates.

Arguments:
- longitude: 2D array of longitude values in radians
- latitude: 2D array of latitude values in radians
- time: DateTime object specifying the time
- density: reference density (default=1)

Returns:
- 2D array of tidal potential values
"""
function compute_potential(λ, φ, time; reference_density=1020)
    # Get celestial body positions
    λ_sun,  φ_lat, R_sun  = get_body_position("SUN", time)
    λ_moon, φ_lat, R_moon = get_body_position("MOON", time)

    # Calculate zenith cosines
    μ_sun  = calculate_zenith_cosine(λ, φ, λ_sun, φ_lat)
    μ_moon = calculate_zenith_cosine(λ, φ, λ_moon, φ_lat)

    # Calculate parallaxes
    χ_sun  = EARTH_RADIUS / R_sun
    χ_moon = EARTH_RADIUS / R_moon

    # Compute potentials
    F_sun  = - G_sun  / R_sun  * χ_sun^2  * Pl(μ_sun, 2)
    F_moon = - G_moon / R_moon * χ_moon^2 * Pl(μ_moon, 2)

    # Apply solid earth tide correction
    ϵ = 1 + LOVE_K2 - LOVE_H2 # = 0.69
    ρₒ = reference_density

    # Return combined, corrected potential
    return ϵ * ρₒ * (F_sun + F_moon)
end

# Create Oceananigans grid
grid = LatitudeLongitudeGrid(size = (180, 80),
                             longitude = (0, 360),
                             latitude = (-80, 80),
                             topology = (Periodic, Bounded, Flat))

# Convert grid coordinates to radians
λ = λnodes(grid, Center(), Center(), Center())
φ = φnodes(grid, Center(), Center(), Center())

# Create field to store tidal potential
Φ = Field{Center, Center, Nothing}(grid)

# Compute potential at current time
t = DateTime(1993, 1, 1, 12)

Nx, Ny, Nz = size(grid)
for i = 1:Nx, j=1:Ny
    Φi = compute_potential(λ[i], φ[j], t)
    Φ[i, j, 1] = Φi
end


greet() = print("Hello World!")

end # module Tidejinks
