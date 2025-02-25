module Tidejinks

using Scratch
using Downloads
using ClimaOcean.DataWrangling: download_progress

kerneldata = Dict(
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

spice_cache::String = ""

function __init__()
    global spice_cache = @get_scratch!("SPICE")

    for (path, url) in kerneldata
        filepath = joinpath(spice_cache, path)
        @info "Downloading file to $filepath..."
        !isfile(filepath) && Downloads.download(url, filepath; progress=download_progress)
    end
end

#####
##### Here we copy the data in scratch into a local directory.
#####
##### This is messy, but _I think_ required to work around the limiations
##### of how SPICE ingests "kernels"; ie the paths in meta_kernel.txt
##### must be local paths? Something to look into.
#####

using SPICE

function wrangle_spice_kernels(metafile="kernels.txt", kernel_relpath="")
                        
    kernel_relpath != "" && !ispath(kernel_relpath) && mkdir(kernel_relpath)

    for path in keys(kerneldata)
        scratch_path = joinpath(spice_cache, path)
        local_path = joinpath(kernel_relpath, path)
        !isfile(local_path) && cp(scratch_path, kernel_relpath)
    end

    if kernel_relpath != "" && kernel_relpath[end] != '/'
        kernel_relpath *= '/'
    end

    metatxt = """
    \\begindata

    KERNELS_TO_LOAD=(
        '$(kernel_relpath)naif0012.tls',
        '$(kernel_relpath)latest_leapseconds.tls',
        '$(kernel_relpath)earth_assoc_itrf93.tf',
        '$(kernel_relpath)de440.bsp',
        '$(kernel_relpath)pck00010.tpc',
        '$(kernel_relpath)gm_de431.tpc',
        '$(kernel_relpath)earth_000101_220503_220207.bpc',
        '$(kernel_relpath)earth_720101_070426.bpc'
    )

    \\begintext
    """

    isfile(metafile) && rm(metafile)
    open(metafile, "w") do file
        write(file, metatxt)
    end

    return nothing
end

# Constants
const EARTH_RADIUS = 6371e3
const LOVE_H2 = 0.61
const LOVE_K2 = 0.30

using Dates

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

"""'h' for 'hack'"""
@inline hsind(φ) = sin(φ * π / 180)
@inline hcosd(φ) = cos(φ * π / 180)

"""
    calculate_zenith_cosine(λ₁, φ₁, λ₂, φ₂)

Calculate cosine of zenith angle using spherical law of cosines.
"""
@inline calculate_zenith_cosine(λ₁, φ₁, λ₂, φ₂) = hsind(φ₁) * hsind(φ₂) + hcosd(φ₁) * hcosd(φ₂) * hcosd(λ₂ - λ₁)

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
function gravitational_parameters()
    # Get gravitational parameters (in m)
    G_sun  = first(bodvrd("SUN", "GM", 1)) * 1e9
    G_moon = first(bodvrd("MOON", "GM", 1)) * 1e9
    return G_sun, G_moon
end

function celestial_positions(time)
    # Get celestial body positions
    X_sun  = get_body_position("SUN", time)
    X_moon = get_body_position("MOON", time)

    return X_sun, X_moon
end

function compute_tidal_potential(λ, φ, time)
    G_sun, G_moon = gravitational_parameters()
    X_sun, X_moon = celestial_positions(time)
    return compute_tidal_potential(λ, φ, X_sun, X_moon, G_sun, G_moon)
end

using LegendrePolynomials

@inline function compute_tidal_potential(λ, φ, X_sun, X_moon, G_sun, G_moon)
    λ_sun,  φ_lat, R_sun  = X_sun
    λ_moon, φ_lat, R_moon = X_moon
    
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

    # Return combined, corrected potential
    return ϵ * (F_sun + F_moon)
end

using Oceananigans
using Oceananigans.Grids: λnode, φnode
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index

const XYField = Field{<:Any, <:Any, Nothing}

function compute_tidal_potential!(Φ::XYField, time)
    G_sun, G_moon = gravitational_parameters()
    X_sun, X_moon = celestial_positions(time)
    parameters = (; G_sun, G_moon, X_sun, X_moon)

    LX, LY, LZ = Oceananigans.Fields.location(Φ)
    ℓx = LX()
    ℓy = LY()
    grid = Φ.grid
    arch = Oceananigans.Architectures.architecture(grid)
    launch!(arch, grid, :xy, _compute_tidal_potential, Φ, grid, ℓx, ℓy, parameters)

    return nothing
end

@kernel function _compute_tidal_potential(Φ, grid, ℓx, ℓy, p)
    i, j = @index(Global, NTuple)
    λ = λnode(i, j, 1, grid, ℓx, ℓy, nothing)
    φ = φnode(i, j, 1, grid, ℓx, ℓy, nothing)
    @inbounds Φ[i, j, 1] = compute_tidal_potential(λ, φ,
                                                   p.X_sun, p.X_moon,
                                                   p.G_sun, p.G_moon)
end

end # module Tidejinks
