# Utility functions, mostly for unit conversions.
# NMR quantities in Bruker notation [1]:
# BF1: Basic transmitter frequency (MHz)
# SF: Spectrometer frequency (MHz)
# SR: BrukerSpectrum reference frequency (Hz)
# SR = (SF - BF1) * 10^6
# SFO1: Transmitter frequency (MHz)
# O1: Transmitter frequency offset (Hz)
# O1 = (SFO1 - BF1) * 10^6
# Ω: Absolute frequency (MHz)
# ω: Frequency relative to SF1 (Hz)
# δ: Chemical shift δ (ppm)
# δ = 10^6 * (Ω - SF) / SF = ω / SF
# 1. Numerical indices indicate different nuclei


limits(s :: S) where S <: AbstractSpectrum = limits(s["O1P"], s["SW"], s["SF"], s["BF1"])
function limits(o1p, sw, sf, bf) 
    shift = 1e6(bf-sf)/bf
    (o1p - sw/2 + shift, o1p + sw/2 + shift)
end

function in(δ, s::S) where S <: AbstractSpectrum
    lo, hi = limits(s)
    lo < δ < hi
end

function chemical_shifts(s::S) where S <: AbstractSpectrum 
    lo, hi = limits(s)
    range(hi; stop=lo, length=length(s[:]))
end


"""
    ppmtomhz_abs(δ, bf[, sr])
Convert chemical shift δ to absolute frequency Ω in MHz.
bf in MHz and sr in Hz.
"""
function ppmtomhz_abs(δ, bf, sr = 0.0)
    bf = bf * 1e6
    1e-6(bf + sr)δ + bf + sr
end

"""
    ppmtohz(δ, bf[, sr])
Convert chemical shift δ in ppm to relative frequency ω in Hz.
bf in MHz and sr in Hz.
"""
function ppmtohz(δ, bf, sr = 0.0)
    (bf + sr*1e-6)δ
end


""" 
    hztoppm(ω, bf[, sr])
Convert relative frequency ω in Hz to chemical shift δ in ppm.
bf in MHz and sr in Hz.
"""
function hztoppm(ω, bf, sr = 0.0)
    ω / (bf + sr*1e-6)
end

"""
    sr(sf, bf)
Return spectral reference (SR in Bruker notation) in Hz.
sf and bf in MHz.
"""
sr(sf, bf) = 1e6(sf - bf)

"""
    ppmtoindex(<:AbstractSpectrum, δ)
Return the index in the processed spectrum corresponding to
chemical shift δ.
"""
function ppmtoindex(s :: S, δ) where S <: AbstractSpectrum
    min_δ, max_δ = limits(s)
    @. Int(cld(s["SI"]*(max_δ - δ), (max_δ - min_δ)))
end
function ppmtoindex(w :: W, δ) where W <: FracSpectrum
    min_δ, max_δ = limits(w)
    Int(cld(length(w.fs)*(max_δ - δ), (max_δ - min_δ)))
end

"""
    ppmtoindex(::PoptSpectrum, δ)
Return the index in the POPT spectrum corresponding to
chemical shift δ. In POPT, this is based on TD, not SI.
"""
function ppmtoindex(p :: P, δ) where P <: PoptSpectrum
    min_δ, max_δ = limits(p)
    @. Int(cld(p["TD"]*(max_δ - δ), (max_δ - min_δ)))
end

function ppmtoindex(s :: S, rng::Tuple{Float64,Float64}) where S <: AbstractSpectrum
    r1, r2 = rng
    if r1>r2
        ppmtoindex(s, r1):ppmtoindex(s, r2)
    else
        ppmtoindex(s, r2):ppmtoindex(s, r1)
    end
end

hztoindex(f, sw, sf, si) = Int(cld(f*si, sw*sf))
hztoindex(s :: S, f) where S <: AbstractSpectrum = hztoindex(f, s["SW"], s["SF"], s["SI"])

Base.length(s :: S) where S <: AbstractSpectrum = length(s[:])

function union_range(ss::AbstractArray{BrukerSpectrum})
    lims = [limits(s) for s in ss]
    l = minimum(lim[1] for lim in lims)
    h = maximum(lim[2] for lim in lims)
    h,l
end

function union_shifts(ss::AbstractArray{BrukerSpectrum})
    h,l = union_range(ss)
    res = minimum(freq_resolution(s) for s in ss)
    h:-res:l
end

union_range(s::BrukerSpectrum) = union_range([s])
union_shifts(s::BrukerSpectrum) = union_shifts([s])

freq_resolution(s::BrukerSpectrum) = s["SW_h"] / s["TD"]

title(s::BrukerSpectrum) = s[s.default_proc].title
title(s::PoptSpectrum) = s.proc.title

# Returns a copy of s with all but the given
# ranges zeroed out. The default proc for s will
# have its intrng adjusted to Δs as well.
function extract(s::BrukerSpectrum, Δs::AbstractArray{Intrng})
    res = deepcopy(s)
    res[res.default_proc].re_ft = zeros(length(s))
    res[res.default_proc].im_ft = zeros(length(s))
    res[res.default_proc].intrng = Δs
    for Δ in Δs
        rng = ppmtoindex(s, Δ)
        res[res.default_proc].re_ft[rng] .= s[s.default_proc].re_ft[rng]
        res[res.default_proc].im_ft[rng] .= s[s.default_proc].im_ft[rng]
    end
    res
end

extract(s::BrukerSpectrum, Δ::Intrng) = extract(s, [Δ])

"""Tunes `param` until `expr` evaluates to zero within δ.
`expr` must be monotonically increasing in terms of `param`."""
macro binary_opt(expr, param, min, max, δ, nattempts=5)
    :(m = $(esc(min));M = $(esc(max));d=$(esc(δ));z=zero(d);n=$(esc(nattempts));
    for i=1:n
        $(esc(param)) = (m+M)/2
        e = $(esc(expr))
        # println("val: $e, param:$((m+M)/2), min: $m, max: $M")
        if norm(e) < d
            break
        elseif e < z
            m = $(esc(param))
        else
            M = $(esc(param))
        end
    end)
end

### Pulse power profile

"""     powerprofile(pulse, dt, sfo1)
----------------------------------------------------------------------------
Pulse power profile for arbitrary pulse shape.
pulse: Amplitude profile of pulse.
dt: 'Sampling frequency' — unit of time between points in time. (s)
sfo1: Centre frequency of pulse. (Hz)"""
function powerprofile(pulse, dt, sfo1)
    ft = abs2.(fft(pulse))
    l = length(pulse)
    df = 1/(l*dt)
    fn = interpolate(ft[1:div(l,2)], BSpline(Cubic(Natural())), OnGrid())
    max = fn[1]
    f -> fn[abs(f-sfo1)/df+1]/max
end


"""    hztoZe(Δω, g)
----------------------------------------------------------------------------
Compute the distance in cm between two frequency points in slice-selective 
experiments, assuming applied gradient g' = g⋅Gmax is accurate / uniform.
In practise the applied gradient is nonuniform, and so this yields distances 
in a mixed "space-frequency" variable ζ.
"""
hztoζ(Δω, g) = Δω / (1e4 * g) 

# some functions to extract the x-axis from a spectrum
ppmaxis(s::S; SI = s["SI"]) where S <: AbstractSpectrum = LinRange(limits(s)..., SI)
hzaxis(s::S; SI = s["SI"]) where S <: AbstractSpectrum = LinRange(s["O1"] - s["SW_h"] / 2, s["O1"] + s["SW_h"] / 2, SI)
ζaxis(s::S; SI = s["SI"]) where S <: AbstractSpectrum = hztoζ.(hzaxis(s; SI = SI))

ppmaxis(w::W; SI = length(w.fs)) where W <: FracSpectrum = LinRange(limits(w)..., SI)
hzaxis(w::W; SI = length(w.fs)) where W <: FracSpectrum = LinRange(w["O1"] - w["SW_h"] / 2, w["O1"] + w["SW_h"] / 2, SI)
ζaxis(w::W; SI = length(w.fs)) where W <: FracSpectrum = hztoζ.(hzaxis(w; SI = SI))

#! TODO implement cmaxis. This function will require a probe-specific gradient map!
# 
# z(ζ) = (2π/γ) × ∫ 1/g(ζ) dζ

cmaxis(s::S; SI = s["SI"]) where S <: AbstractSpectrum = @error "cmaxis() not implemented yet."

"""
    trim_spectrum(S, δ) -> (idx1, idx2)

Ends of the spectrum are often distorted by DSP, and look ugly in plots or have weird phase errors.
-   δ   a peak in ppm, if given will indices for the region δ +- 2 Hz 
"""
function trim_spectrum(s :: S, δ = (); pc = 0.01) where S <: AbstractSpectrum
    # default trim, short back & sides
    if isempty(δ)
        n = size(s)[1]
        idx1, idx2 = ceil(Int, pc * n), floor(Int, (1 - pc)n)
     
    # Otherwise phase around a specific peak.
    else 
        idx1 = ppmtoindex(s, δ + 3)
        idx2 = ppmtoindex(s, δ - 3)
    end
    return (idx1, idx2)
end
export ppmaxis, hzaxis, ppmtomhz_abs, ppmtohz, ppmtoindex, hztoindex, hztoppm, hztoζ, union_range, union_shifts, limits, trim_spectrum
