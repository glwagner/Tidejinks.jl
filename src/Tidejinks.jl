module Tidejinks

using Scratch
using Downloads
using LegendrePolynomials
using Oceananigans
using ClimaOcean.DataWrangling: download_progress

# URLs to NAIF data
# See https://naif.jpl.nasa.gov/pub/naif for more information.
const NAIF = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels"

kerneldata = Dict(
    # "Leap seconds kernel
    "de440.bsp"                      => NAIF * "/spk/planets/de440.bsp",
    "earth_assoc_itrf93.tf"          => NAIF * "/fk/planets/earth_assoc_itrf93.tf",
    "earth_000101_220503_220207.bpc" => NAIF * "/pck/earth_000101_220503_220207.bpc",
    "earth_720101_070426.bpc"        => NAIF * "/pck/earth_720101_070426.bpc",
    "earth_000101_220503_220207.bpc" => NAIF * "/pck/earth_000101_220503_220207.bpc",
    "gm_de431.tpc"                   => NAIF * "/pck/gm_de431.bpc",
    "pck00010.tpc"                   => NAIF * "/pck/pck00010.tpc",
    "latest_leapseconds.tls"         => NAIF * "/lsk/latest_leapseconds.tls",
    # "Leap seconds kernel
    "naif0012.tls"                   => NAIF * "/lsk/naif0012.tls",
)

backup_kerneldata = Dict(
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

"""
    retrying_download(url::String, dest::String; retries=3, delay=2)

Attempts to download a file from `url` to `dest`.
Retries up to `max_retries` times with an increasing delay.
Returns `true` if successful, `false` otherwise.
"""
function retrying_download(url::String, dest::String; progress=nothing, retries=3, delay=4)
    # Downloads.download(url, dest; progress)

    if !isfile(dest)
        for attempt in 1:retries
            try
                Downloads.download(url, dest; progress)
                @info "Successfully downloaded $url -> $dest"
                return true
            catch e
                @warn "Download failed (attempt $attempt/$retries) for $url: $e"
                if attempt < retries
                    sleep(delay * attempt) # Exponential backoff
                    @warn "Retrying in $(delay * attempt) seconds..."
                end
            end
        end

        @warn "Failed to download $url after $retries attempts"
        return false
    else
        return true
    end
end

function __init__()
    global spice_cache = @get_scratch!("SPICE")

    for (path, url) in kerneldata
        filepath = joinpath(spice_cache, path)
        @info "Downloading file to $filepath..."
        success = retrying_download(url, filepath; progress=download_progress)

        if !success # we have a backup!
            backup_url = backup_kerneldata[path]
            @info "Download failed for $url."
            @info "Trying backup_url $backup_url..."
            retrying_download(backup_url, filepath; progress=download_progress)
        end
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
        !isfile(local_path) && cp(scratch_path, local_path)
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

"""
    calculate_zenith_cosine(λ₁, φ₁, λ₂, φ₂)

Calculate cosine of zenith angle using spherical law of cosines.
"""
@inline calculate_zenith_cosine(λ₁, φ₁, λ₂, φ₂) = sin(φ₁) * sin(φ₂) + cos(φ₁) * cos(φ₂) * cos(λ₂ - λ₁)

"""
    gravitational_parameters()

Return `G_sun, G_moon`, the gravitational constants for the sun and moon,
respectively.

Uses the CSPICE function `bodvrd`.
See https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html.
"""
function gravitational_parameters()
    # Get gravitational parameters (in m)
    # https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html
    G_sun  = first(bodvrd("SUN", "GM", 1)) * 1e9
    G_moon = first(bodvrd("MOON", "GM", 1)) * 1e9
    return G_sun, G_moon
end

"""
    celestrial_positions(time)

Return `X_sun, X_moon`, 3-tuples that represent the
longitude, latitude, and distance of the sun and moon
respectively relative to Earth.
"""
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

"""
    compute_tidal_potential(λ, φ, X_sun, X_moon, G_sun, G_moon)

Compute the effective tidal potential corrected for the
solid Earth tide at Earth longitude `λ` and Earth latitude `φ`,
given spherical position of the sun and moon `X_sun` and `X_moon`,
and the gravitational constants for the sun and moon `G_sun` and `G_moon`.
"""
@inline function compute_tidal_potential(λ, φ, X_sun, X_moon, G_sun, G_moon)
    λ_sun,  φ_lat, R_sun  = X_sun
    λ_moon, φ_lat, R_moon = X_moon
    
    # Calculate zenith cosines
    μ_sun  = calculate_zenith_cosine(λ, φ, λ_sun, φ_lat)
    μ_moon = calculate_zenith_cosine(λ, φ, λ_moon, φ_lat)

    !(-1 < μ_sun < 1)  && @warn("μ_sun=$μ_sun lies outside [-1, 1]")
    !(-1 < μ_moon < 1) && @warn("μ_moon=$μ_moon lies outside [-1, 1]")

    μ_sun  = clamp(μ_sun,  -1, 1)
    μ_moon = clamp(μ_moon, -1, 1)

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

"""
    compute_tidal_potential!(Φ, time)

Compute the effective tidal potential corrected for the
solid Earth tide at `time` and over the nodes of `Φ`.
"""
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

    # Convert λ to -180..180
    λ -= 360 * (λ > 180)

    # Convert to radians
    λ *= π/180 
    φ *= π/180 

    @inbounds Φ[i, j, 1] = compute_tidal_potential(λ, φ,
                                                   p.X_sun, p.X_moon,
                                                   p.G_sun, p.G_moon)
end

end # module Tidejinks

