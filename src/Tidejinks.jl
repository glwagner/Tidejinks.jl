module Tidejinks

using Scratch
using Downloads
using LegendrePolynomials
using Oceananigans

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
    "gm_de431.tpc"                   => NAIF * "/pck/gm_de431.tpc",
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
    retrying_download(url::String, dest::String; max_retries=3, delay=2.0)

Attempts to download a file from `url` to `dest`.
Retries up to `max_retries` times with an increasing delay.
Returns `true` if successful, `false` otherwise.
"""
function retrying_download(url::String, dest::String; max_retries=3, delay=2)
    if !isfile(dest)
        for attempt in 1:max_retries
            try
                Downloads.download(url, dest)
                @info "Successfully downloaded $url -> $dest"
                return true
            catch e
                @warn "Download failed (attempt $attempt/$max_retries) for $url: $e"
                if attempt < max_retries
                    sleep(delay * attempt)  # Exponential backoff
                    @warn "Retrying in $(delay * attempt) seconds..."
                end
            end
        end

        @warn "Failed to download $url after $max_retries attempts"
    end

    return false
end

function __init__()
    global spice_cache = @get_scratch!("SPICE")

    for (path, url) in kerneldata
        filepath = joinpath(spice_cache, path)
        @info "Downloading file to $filepath..."
        success = retrying_download(url, filepath)
        if !success
            backup_url = backup_kerneldata[path]
            @info "Download failed for $url."
            @info "Trying backup_url $backup_url..."
            retrying_download(backup_url, filepath)
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
    get_body_position(body, time, reference_frame="ITRF93")

Get longitude, latitude, and distance of celestial body in the Earth-fixed `reference_frame`.
"""
function get_body_position(body, time, reference_frame="ITRF93")
    time_str = Dates.format(time, "yyyy-mm-dd HH:MM:SS") * " UTC"
    ephemeris_time = str2et(time_str)
    position_rectangular, _ = spkpos(body, ephemeris_time, reference_frame, "NONE", "EARTH")
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
function celestial_positions(time, reference_frame="ITRF93")
    # Get celestial body positions
    X_sun  = get_body_position("SUN", time, reference_frame)
    X_moon = get_body_position("MOON", time, reference_frame)

    return X_sun, X_moon
end

function compute_tidal_potential(λ, φ, time, reference_frame="ITRF93")
    G_sun, G_moon = gravitational_parameters()
    X_sun, X_moon = celestial_positions(time, reference_frame)
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
    λ_sun,  φ_sun,  R_sun  = X_sun
    λ_moon, φ_moon, R_moon = X_moon

    # Calculate zenith cosines
    μ_sun  = calculate_zenith_cosine(λ, φ, λ_sun, φ_sun)
    μ_moon = calculate_zenith_cosine(λ, φ, λ_moon, φ_moon)

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
    compute_tidal_potential!(Φ, time, reference_frame="ITRF93")

Compute the effective tidal potential corrected for the
solid Earth tide at `time` and over the nodes of `Φ`.
"""
function compute_tidal_potential!(Φ::XYField, time, reference_frame="ITRF93")
    G_sun, G_moon = gravitational_parameters()
    X_sun, X_moon = celestial_positions(time, reference_frame)
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

#####
##### The tidal potential as a function of model time, for evaluation in kernels
#####

using Adapt
using Oceananigans.Architectures: architecture, on_architecture

# Rotation rate of the IAU_EARTH model, 360.9856235ᵒ per day
const EARTH_ROTATION_RATE = 2π * 360.9856235 / 360 / 86400

struct TidalPotential{FT, C, P} <: Function
    first_sample_time :: FT
    sampling_interval :: FT
    rotation_rate :: FT
    coefficients :: C
    phases :: P
end

Adapt.@adapt_structure TidalPotential

"""
    TidalPotential(grid, dates; epoch = first(dates), sampling_interval = Hour(1), reference_frame = "ITRF93")

Return the tidal potential `Φ(λ, φ, t)` (m² s⁻²) over `dates = (start_date, end_date)`,
where `λ` and `φ` are longitude and latitude in degrees and `t` is the model time
in seconds since `epoch`. `Φ` is the potential of [`compute_tidal_potential`](@ref),
written with the addition theorem as

    Φ = D₀ P₂(sin φ) + sin φ cos φ (D₁ᶜ cos χ + D₁ˢ sin χ) + cos²φ (D₂ᶜ cos 2χ + D₂ˢ sin 2χ),

where `χ = λ + Ω t` and `Ω` is Earth's rotation rate. The five coefficients `D` vary on the
time scales of the lunar and solar orbits. They are computed from the positions of the Sun
and the Moon every `sampling_interval` and interpolated in time with four-point Lagrange
polynomials. The rotation phase `Ω t` is stored modulo 2π at each sample, so `Φ` stays
accurate when `t` is single precision.

SPICE kernels with the ephemerides and the `reference_frame` must be loaded.
"""
function TidalPotential(grid, dates; epoch = first(dates), sampling_interval = Hour(1), reference_frame = "ITRF93")
    start_date, end_date = dates
    G_sun, G_moon = gravitational_parameters()
    Δτ = Millisecond(sampling_interval)
    sample_dates = (start_date - Δτ):Δτ:(end_date + 2Δτ)
    model_time(date) = Dates.value(date - epoch) / 1000

    coefficients = zeros(5, length(sample_dates))
    phases = zeros(length(sample_dates))
    for (n, date) in enumerate(sample_dates)
        ϑ = EARTH_ROTATION_RATE * model_time(date)
        coefficients[:, n] .= rotate(tidal_coefficients(date, G_sun, G_moon, reference_frame), ϑ)
        phases[n] = mod(ϑ, 2π)
    end

    FT = eltype(grid)
    arch = architecture(grid)
    return TidalPotential(convert(FT, model_time(first(sample_dates))),
                          convert(FT, Dates.value(Δτ) / 1000),
                          convert(FT, EARTH_ROTATION_RATE),
                          on_architecture(arch, convert(Array{FT}, coefficients)),
                          on_architecture(arch, convert(Array{FT}, phases)))
end

# Coefficients of P₂(sin φ), sin φ cos φ (cos λ, sin λ), and cos²φ (cos 2λ, sin 2λ) in the potential at `date`
function tidal_coefficients(date, G_sun, G_moon, reference_frame)
    ϵ = 1 + LOVE_K2 - LOVE_H2
    C = zeros(5)
    for (G, (λ, φ, R)) in zip((G_sun, G_moon), celestial_positions(date, reference_frame))
        A = - ϵ * G / R * (EARTH_RADIUS / R)^2
        s, c = sin(φ), cos(φ)
        C .+= A .* ((3s^2 - 1) / 2, 3s * c * cos(λ), 3s * c * sin(λ), 3c^2 / 4 * cos(2λ), 3c^2 / 4 * sin(2λ))
    end
    return C
end

# Coefficients of the same potential in terms of χ = λ + ϑ
function rotate(C, ϑ)
    c₁, s₁ = cos(ϑ), sin(ϑ)
    c₂, s₂ = cos(2ϑ), sin(2ϑ)
    return (C[1], C[2] * c₁ - C[3] * s₁, C[3] * c₁ + C[2] * s₁, C[4] * c₂ - C[5] * s₂, C[5] * c₂ + C[4] * s₂)
end

# Four-point Lagrange interpolation at θ ∈ [0, 1] between y₀ and y₁
@inline lagrange(y₋₁, y₀, y₁, y₂, θ) = - θ * (θ - 1) * (θ - 2) / 6 * y₋₁ + (θ + 1) * (θ - 1) * (θ - 2) / 2 * y₀ -
                                         (θ + 1) * θ * (θ - 2) / 2 * y₁ + (θ + 1) * θ * (θ - 1) / 6 * y₂

@inline function (Φ::TidalPotential)(λ, φ, t)
    Δτ = Φ.sampling_interval
    s = (t - Φ.first_sample_time) / Δτ
    N = length(Φ.phases)
    n = clamp(unsafe_trunc(Int, floor(s)) + 1, 2, N - 2)
    θ = s - (n - 1)
    D = Φ.coefficients

    @inbounds begin
        D₀  = lagrange(D[1, n-1], D[1, n], D[1, n+1], D[1, n+2], θ)
        D₁ᶜ = lagrange(D[2, n-1], D[2, n], D[2, n+1], D[2, n+2], θ)
        D₁ˢ = lagrange(D[3, n-1], D[3, n], D[3, n+1], D[3, n+2], θ)
        D₂ᶜ = lagrange(D[4, n-1], D[4, n], D[4, n+1], D[4, n+2], θ)
        D₂ˢ = lagrange(D[5, n-1], D[5, n], D[5, n+1], D[5, n+2], θ)
        χ = deg2rad(λ) + Φ.phases[n] + Φ.rotation_rate * (t - Φ.first_sample_time - (n - 1) * Δτ)
    end

    μ = sin(deg2rad(φ))
    c = cos(deg2rad(φ))

    return D₀ * (3μ^2 - 1) / 2 + μ * c * (D₁ᶜ * cos(χ) + D₁ˢ * sin(χ)) + c^2 * (D₂ᶜ * cos(2χ) + D₂ˢ * sin(2χ))
end

end # module Tidejinks

