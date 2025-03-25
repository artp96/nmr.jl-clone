# Some utilities
const S = String
"Integration range in the frequency domain, values given in ppm."
const Intrng{T} = Tuple{T, T} where T <: AbstractFloat
abstract type AbstractSpectrum end

"""
Processed (frequency domain) NMR spectrum. Parameters follow the Bruker conventions:
- **re_ft**: real part of Fourier transform
- **im_ft**: imaginary part of Fourier transform
- **params**: dictionary of processing parameters, mostly from `proc` file for Bruker data
- **intrng**: list of integration ranges
- **title**: Spectrum title
"""

mutable struct ProcessedSpectrum{T <: AbstractFloat, V <: AbstractVector{T}, S <: String} <: AbstractSpectrum
    re_ft :: V
    im_ft :: V
    params :: Dict{S, Any}
    intrng :: Union{Intrng{T},Missing}
    procno :: Int
    title :: S
    function ProcessedSpectrum{T, V, S}(re :: V, im :: V, par :: Dict{S, Any}, intrng :: Union{Intrng{T}, Missing}, pn :: Int, t :: S) where {
        T <: AbstractFloat, V <: AbstractVector{T}, S <: String
        }
        # Validate that `intrng` is a Tuple of two elements
        
        ismissing(intrng) && return new{T, V, S}(re, im, par, missing, pn, t) 

        if !(length(intrng) == 2)
            throw(ArgumentError("`intrng` must be a tuple of length 2"))
        end
        return new{T, typeof(re), S}(re, im, par, intrng, pn, t)
    end
end

 """
    Outer Constructor Definition"""
function ProcessedSpectrum(re::V, im::V, par::Dict{S, Any}, intrng::Union{Intrng{T}, Missing}, pn::Int, t::S) where {
    T <: AbstractFloat, V <: AbstractVector{T}, S <: String
}
    ProcessedSpectrum{T, V, S}(re, im, par, intrng, pn, t)
end

Base.getindex(p::ProcessedSpectrum, n::Int) = p.intrng[n]
Base.getindex(p::ProcessedSpectrum, param::AbstractString) = p.params[param]
Base.getindex(p::ProcessedSpectrum, ::Colon) = p.re_ft
Base.getindex(p::ProcessedSpectrum, a::AbstractArray) = p.re_ft[a]
Base.view(p::ProcessedSpectrum, v) = @view p.re_ft[v]

"""
Unprocessed NMR experiment: FID and any associated processed spectra.

Following Bruker convention:
- **fid**: Acquired free induction decay (FID) signal
- **acqu**: Acquisition parameters, mostly from `acqu` file for Bruker data
- **procs**: Dictionary of associated processed spectra
- **default_proc**: Default processed spectrum to be used in operations that require a processed spectrum
- **name**: Experiment name
- **expno**: Experiment number if, e.g., part of a Bruker dataset
"""
mutable struct Spectrum{T <: AbstractFloat, V <: AbstractVector{<: Complex{T}}, P <: ProcessedSpectrum{T}}  <: AbstractSpectrum
    fid :: V
    acqu :: Dict{String, Any}
    procs :: Dict{Int, P}
    default_proc :: Int
    name :: String
    expno :: Int
        
    # Inner constructor to handle weakly typed dicts & split the FID.
    function Spectrum(f :: V, a :: Dict{S, Any}, p :: Dict{I, P}, d :: I, n :: String, e :: I) where {
        # Type-fu is intensifying #
        T <: AbstractFloat, 
        V <: AbstractVector{T}, 
        P <: ProcessedSpectrum{T, V}, 
        S <: AbstractString, 
        I <: Int}
        # process the FID, {ℝ, ℝ} -> ℂ
        f = split_fid(f)  # Ensure that the transformed FID is of the correct type
        if !(f isa AbstractVector{<: Complex})
            throw(ArgumentError("The split FID must be of type AbstractVector{$C}, but got $(typeof(f))"))
        end
        W = typeof(f)
        return new{T, W, P}(f, a, p, d, n, e)
    end
end

Base.getindex(s::Spectrum, i::Int) = s.procs[i]
Base.getindex(s::Spectrum, ::Colon) = s.procs[s.default_proc].re_ft
Base.getindex(s::Spectrum, a::AbstractArray) = s.procs[s.default_proc].re_ft[a]
Base.getindex(s::Spectrum, rng::Tuple{Float64,Float64}) = s[ppmtoindex(s,rng)]
Base.getindex(s::Spectrum, δ::Float64) = s[s.default_proc].re_ft[ppmtoindex(s, δ)]
Base.view(s::Spectrum, v) = @view s.procs[s.default_proc][v]

function Base.getindex(s::Spectrum, param::AbstractString)
    try
        s.acqu[param]
    catch err
        s[s.default_proc][param]
    end
end
Base.setindex!(s::Spectrum, d::AbstractArray, ::Colon) = (s[s.default_proc].re_ft .= d)
Base.setindex!(s::Spectrum, d::AbstractArray, r::UnitRange) = (s[s.default_proc].re_ft[r] .= d)
Base.setindex!(s::Spectrum, d, rng::Tuple{Float64, Float64}) = (s[s.default_proc].re_ft[ppmtoindex(s,rng)]=d)

function Spectrum(fid :: V, acqu :: Dict{S, Any}, proc :: P) where {
    T <: AbstractFloat, V <: AbstractVector{T}, S <: AbstractString, P <: ProcessedSpectrum{T, V}
    }
    return Spectrum(fid, acqu, Dict{Int, ProcessedSpectrum{T, V}}(1=>proc), 1, "", "")
end

"""
    PoptSpectrum(f :: V, a :: Dict{S, Any}, p :: Int, e:: Int,  n :: S, fp :: S, prot :: DF, ser :: V, proc :: P) where {
        T <: AbstractFloat,
        C <: Complex{T},
        V <: AbstractVector{C},
        A <: AbstractArray{C}
        }

POPT array of 1-D spectra.

Following Bruker convention:
- **fid**: Complex FID from the top-level expno 
- **acqu**: Acquisition parameters, mostly from `acqu` file for Bruker data
- **name**: Experiment name
- **expno**: Experiment number if, e.g., part of a Bruker dataset
- **fpath**: File path should point to the correct protocol with numerical suffix
- **protocol**: Processed POPT protocol file
- **ser**: N-dimensional full spectrum from saved ser, may not always exist 
- **proc** Bruker-processed POPT spectrum (f1p ≤ δ ≤ f2p).
- **vars**: Dict{S, Vector{T}} - array of POPT variables

 ***TODO*** 
 1. get working for N-dimensional spectra
 2. gpu-compatible with flux.@functor ?
"""
mutable struct PoptSpectrum{
        T <: AbstractFloat,
        C <: Complex{T},
        V <: AbstractVector{C},
        A <: AbstractArray{C}
        } <: AbstractSpectrum

    fid :: V # should probably shadow the ser if it exists
    acqu :: Dict{S, Any}
    procno :: Int # for popt, usually 899, 898, etc.
    expno :: Int
    name :: S
    fpath :: S # Should point to the correct protocol with numerical suffix
    protocol :: D where D <: AbstractDataFrame
    ser :: A # raw 2D spectrum, may not always exist
    proc :: ProcessedSpectrum # Bruker processed POPT output, as array
    vars :: Dict{S, Vector{T}}

    """
        PoptSpectrum(f :: V, a :: Dict{S, Any}, p :: Int, e:: Int,  n :: S, fp :: S, prot :: DF, ser :: V, proc :: P) where { 
        T <: AbstractFloat,
        C <: Complex{T},
        V <: AbstractVector{C},
        S <: String,
        DF <: AbstractDataFrame,
        P <: ProcessedSpectrum
        }

Inner constructor to restructure the popt array automatically, & store the popt vars.
Outer constructors should be used to get all the variables from a popt file.
    """
    function PoptSpectrum(f :: V, a :: Dict{S, Any}, p :: Int, e:: Int,  n :: S, fp :: S, prot :: DF, ser :: V, proc :: P) where {
        T <: AbstractFloat,
        C <: Complex{T},
        V <: AbstractVector{C},
        S <: String,
        DF <: AbstractDataFrame,
        P <: ProcessedSpectrum
        }
        dims, vars = get_ArrayPoptDims(prot)
        # Validate if the data length matches expected size
        expected_size = prod(dims)
        actual_size = length(ser)
        length(ser) % prod(dims) != 0 && throw(ArgumentError("
            \n Mismatch between ser length ($actual_size) \n and expected size ($expected_size)."))
        
        # Ensure the array correctly includes each fid as the last dimension;
        # n.b. ÷() is integer divisor operator, /() returns f64.
        # procdims = tuple(dims..., length(proc) ÷ prod(dims))
        # proc = reshape(proc, procdims)
        # Handle the case of no serfile.
        serdims = tuple(length(ser) ÷ prod(dims), dims...)
        ser = isempty(ser) ? Complex(Float64[]) : reshape(ser, serdims)
        A = typeof(ser)
        return new{T, C, V, A}(f, a, p, e, n, fp, prot, ser, proc, vars)
    end

end

"""
    PoptSpectrum(f :: V, a :: Dict{S, Any}, p :: Int, e:: Int,  n :: S, fp :: S, prot :: DF, ser :: V, proc :: P) where {
        T <: AbstractFloat, 
        V <: AbstractVector{T},
        S <: String, 
        DF <: AbstractDataFrame, 
        P <: ProcessedSpectrum }

Outer constructor to convert raw fid/ser to complex fid/ser.
"""
function PoptSpectrum(f :: V, a :: Dict{S, Any}, p :: Int, e:: Int,  n :: S, fp :: S, prot :: DF, ser :: V, proc :: P) where {
        T <: AbstractFloat, V <: AbstractVector{T}, S <: String, DF <: AbstractDataFrame, P <: ProcessedSpectrum
        }
        f = split_fid(f)
        ser = split_fid(ser)
        return PoptSpectrum(f, a, p, e, n, fp, prot, ser, proc)
end

Base.getindex(s::PoptSpectrum, n::Int) = n == 1 ? s.proc : @error("getindex(s, n) not defined for n ≠ 1 on POPT spectra.")

Base.getindex(s::PoptSpectrum, ::Colon) = real(s.ser)
Base.getindex(s::PoptSpectrum, dims::NTuple{Union{Colon,Int}}) = real(s.ser[dims...])
Base.getindex(s::PoptSpectrum, a::AbstractArray) = s.proc[a]
Base.getindex(s::PoptSpectrum, rng::Tuple{Float64,Float64}) = s[ppmtoindex(s,rng)]
# Base.getindex(s::PoptSpectrum, δ::Float64) = s[s.default_proc].re_ft[ppmtoindex(s, δ)]
# Base.view(s::PoptSpectrum, v) = @view s.procs[s.default_proc][v]
function Base.getindex(s::PoptSpectrum, param::AbstractString)
    try
        s.acqu[param]
    catch err
        getindex(s.proc, param)
    end
end


import Base.im
""" 
    synthesize_ser(s :: S) where S <: ProcessedSpectrum

Function to produce a synthetic serfile from POPT-processed data - dubious! 
No trivial way to find units!
May not act correctly if the procno is not a popt procno!
"""
synthesize_ser(s :: S) where S <: ProcessedSpectrum = s.re_ft + base.im * s.im_ft

"""
    split_fid(fid :: V) where {V <: AbstractVector}

Split a fresh FID into it's real and complex parts. Zero-fills to 2N.
N-D spectra are still acquired time-domain sequentially, so this should work for all N-D spectra?
"""
function split_fid(fid :: V) where {V <: AbstractVector}
    fid = complex(fid)
    fid[1:2:end] += 1im * fid[2:2:end]
    # zerofilling
    fid[2:2:end] .= 0 + 0im
    return fid
end

import Base.Threads
"""
    split_fid(ser :: A) where {T, N, A <: AbstractArray{T, N}} 

Split an array of fids, vector-wise.
"""
function split_fid(ser :: A) where {T, N, A <: AbstractArray{T, N}} 
    Threads.@threads for f in eachfibre(ser)
        f = split_fid(f)
    end
    return ser
end

export Spectrum, PoptSpectrum, ProcessedSpectrum, AbstractSpectrum
