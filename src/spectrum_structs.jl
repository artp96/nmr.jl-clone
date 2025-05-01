# Some utilities
const S = String
"Integration range in the frequency domain, values given in ppm."
const Intrng{T} = Tuple{T, T} where T <: AbstractFloat

import Base.complex, Base.imag, Base.real

"""
    Abstract Supertype of spectrum structs.
"""
abstract type AbstractSpectrum end


""" ProcessedSpectrum <: AbstractSpectrum
    - **re_ft**: real part of Fourier transform
    - **im_ft**: imaginary part of Fourier transform
    - **params**: dictionary of processing parameters, mostly from `proc` file for Bruker data
    - **intrng**: list of integration ranges
    - **title**: BrukerSpectrum title
----------------------------------------------------------------------------------------------
Processed (Bruker frequency domain) NMR spectrum. Parameters follow the Bruker conventions:
"""
mutable struct ProcessedSpectrum{T <: AbstractFloat, V <: AbstractVector{T}, S <: String} <: AbstractSpectrum
    re_ft :: V
    im_ft :: V
    params :: Dict{S, Any}
    intrng :: Union{Missing, Vector{Tuple{T, T}}}
    procno :: Int
    title :: S
    function ProcessedSpectrum{T, V, S}(re :: V, im :: V, par :: Dict{S, Any}, intrng :: I, pn :: Int, t :: S) where {
        T <: AbstractFloat, V <: AbstractVector{T}, S <: String, I <: Union{Missing, Vector{Tuple{T, T}} }

        }
        # Validate that `intrng` is a Tuple of two elements
        
        if ismissing(intrng) 
            intrng = Vector{Tuple{T, T}}()
        end
        
        return new{T, V, S}(re, im, par, intrng, pn, t)
    end
end

 """
    ProcessedSpectrum(re::V, im::V, par::D, intrng::I, pn::Int, t::S) where {
    T <: AbstractFloat, 
    V <: AbstractVector{T}, 
    S <: String, 
    I <: Union{ Tuple{T, T}, Vector{Tuple{T, T}}, Missing},
    D <: AbstractDict
    }
 -----------------------------------------------------------------------------------------
    Outer Constructor Definition

    Creates a `ProcessedSpectrum` object, accepting `Intrng{T}` or `Vector{Intrng{T}}`.
 """
function ProcessedSpectrum(re::V, im::V, par::D, intrng::I, pn::Int, t::S) where {
    T <: AbstractFloat, 
    V <: AbstractVector{T}, 
    S <: String, 
    I <: Union{ Tuple{T, T}, Vector{Tuple{T, T}}, Missing},
    D <: AbstractDict
}
    intrng isa Tuple{T, T} ? intrng = [intrng] : nothing
        
    ProcessedSpectrum{T, V, S}(re, im, par, intrng, pn, t)
end

Base.getindex(p::ProcessedSpectrum, n::Int) = p.intrng[n]
Base.getindex(p::ProcessedSpectrum, param::AbstractString) = p.params[param]
Base.getindex(p::ProcessedSpectrum, ::Colon) = p.re_ft
Base.getindex(p::ProcessedSpectrum, a::AbstractArray) = p.re_ft[a]
Base.view(p::ProcessedSpectrum, v) = @view p.re_ft[v]
Base.abs(p::ProcessedSpectrum) = abs.(complex(p))

# Functional shorthand into procno spectra
# may need correction for baseline in every case for Procnos?
Base.real(p::P) where P<:ProcessedSpectrum = p.re_ft #./ p.re_ft[1]
Base.imag(p::P) where P<:ProcessedSpectrum = im * p.im_ft #./ p.im_ft[1]
Base.complex(p::P) where P<:ProcessedSpectrum = p.re_ft .+ im * p.im_ft
Base.size(p :: P) where P <: ProcessedSpectrum = size(p.re_ft)

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
mutable struct BrukerSpectrum{T <: AbstractFloat, V <: AbstractVector{<: Complex{T}}, P <: ProcessedSpectrum{T}}  <: AbstractSpectrum
    fid :: V
    acqu :: Dict{String, Any}
    procs :: Dict{Int, P}
    default_proc :: Int
    name :: String
    expno :: Int
        
    # Inner constructor to handle weakly typed dicts & split the FID.
    function BrukerSpectrum(f :: V, a :: Dict{S, Any}, p :: Dict{I, P}, d :: I, n :: String, e :: I) where {
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
        # Copy the acqupars into each procno
        s = new{T, W, P}(f, a, p, d, n, e)
        ParamDict(s)
        return s
    end
end

Base.getindex(s::BrukerSpectrum, i::Int) = s.procs[i]
Base.getindex(s::BrukerSpectrum, ::Colon) = s.procs[s.default_proc].re_ft
Base.getindex(s::BrukerSpectrum, a::AbstractArray) = s.procs[s.default_proc].re_ft[a]
Base.getindex(s::BrukerSpectrum, rng::Tuple{Float64,Float64}) = s[ppmtoindex(s,rng)]
Base.getindex(s::BrukerSpectrum, δ::Float64) = s[s.default_proc].re_ft[ppmtoindex(s, δ)]
Base.view(s::BrukerSpectrum, v) = @view s.procs[s.default_proc][v]
Base.size(s::BrukerSpectrum) = size(s[s.default_proc].re_ft)
Base.abs(s::BrukerSpectrum) = abs(s[s.default_proc])

function Base.getindex(s::BrukerSpectrum, param::AbstractString)
    try
        s.acqu[param]
    catch err
        s[s.default_proc][param]
    end
end
Base.setindex!(s::BrukerSpectrum, d::AbstractArray, ::Colon) = (s[s.default_proc].re_ft .= d)
Base.setindex!(s::BrukerSpectrum, d::AbstractArray, r::UnitRange) = (s[s.default_proc].re_ft[r] .= d)
Base.setindex!(s::BrukerSpectrum, d, rng::Tuple{Float64, Float64}) = (s[s.default_proc].re_ft[ppmtoindex(s,rng)]=d)

function BrukerSpectrum(fid :: V, acqu :: Dict{S, Any}, proc :: P) where {
    T <: AbstractFloat, V <: AbstractVector{T}, S <: AbstractString, P <: ProcessedSpectrum{T, V}
    }
    s = BrukerSpectrum(fid, acqu, Dict{Int, ProcessedSpectrum{T, V}}(1=>proc), 1, "", "")
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
- **ft**: an fft_plan for efficient on the fly fft of ser
- **acqu**: Acquisition parameters, mostly from `acqu` file for Bruker data
- **name**: Experiment name
- **expno**: Experiment number if, e.g., part of a Bruker dataset
- **fpath**: File path should point to the correct protocol with numerical suffix
- **protocol**: Processed POPT protocol, a named tuple of Popt variables
- **ser**: N-dimensional full spectrum from saved ser, may not always exist 
- **proc** Bruker-processed POPT spectrum, δ ∈ [f1p...f2p].

 ***TODO*** 
 1. get working for N-dimensional spectra
 2. gpu-compatible with flux.@functor ?
"""
mutable struct PoptSpectrum{
        T <: AbstractFloat,
        C <: Complex{T},
        A <: AbstractArray{C}
        } <: AbstractSpectrum

    ft :: FT where FT <: AbstractFFTs.Plan # store for fast fft of this object 
    acqu :: Dict{S, Any}
    procno :: Int # for popt, usually 899, 898, etc.
    expno :: Int
    name :: S
    fpath :: S # Should point to the correct protocol with numerical suffix
    protocol :: NT where NT <: NamedTuple
    ser :: A # transformed 2-D spectrum
    proc :: ProcessedSpectrum # Bruker processed POPT output, as array

    """
        PoptSpectrum(f :: V, a :: Dict{S, Any}, p :: Int, e:: Int,  n :: S, fp :: S, prot :: DF, ser :: V, proc :: P) where { 
        T <: AbstractFloat,
        C <: Complex{T},
        V <: AbstractVector{C},
        S <: String,
        DF <: AbstractDataFrame,
        P <: ProcessedSpectrum
        }

Inner constructor to restructure the popt array automatically, & store the popt vars extracted from a DataFrame.
Outer constructors should be used to get all the variables from a popt file.
    """
    function PoptSpectrum(a :: Dict{S, Any}, p :: Int, e:: Int,  n :: S, fpath :: S, prot :: DF, ser :: V, proc :: P) where {
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
        
        ft = plan_fft(ser, [1])
        
        
        return new{T, C, A}(ft, a, p, e, n, fpath, vars, ser, proc)
    end
end

import FFTW: fft
fft(s :: P) where P <:  PoptSpectrum = s.ft * s.ser

"""
    PoptSpectrum(f :: V, a :: Dict{S, Any}, p :: Int, e:: Int,  n :: S, fp :: S, prot :: DF, ser :: V, proc :: P) where {
        T <: AbstractFloat, 
        V <: AbstractVector{T},
        S <: String, 
        DF <: AbstractDataFrame, 
        P <: ProcessedSpectrum }

Outer constructor to convert raw fid/ser to complex fid/ser.
"""
function PoptSpectrum(a :: Dict{S, Any}, p :: Int, e:: Int,  n :: S, fp :: S, prot :: DF, ser :: V, proc :: P) where {
        T <: AbstractFloat, V <: AbstractVector{T}, S <: String, DF <: AbstractDataFrame, P <: ProcessedSpectrum
        }
        ser = split_fid(ser)
        return PoptSpectrum(a, p, e, n, fp, prot, ser, proc)
end

import Base.size
Base.size(p :: P) where P <: PoptSpectrum = size(p.ser)


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
    return fid[1:2:end]
end

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


"""
    ParamDict(s <: AbstractSpectrum) -> s
- acqupars :: Dict
- procpars :: Dict
----------------------------------------------------------------------------------------
For procnos, stores a copy of the acqupars so these can be referenced directly without
knowledge of the parent procs. Note that this stores the TopSpin parameters, immutably.
Referencing procpars at the BrukerSpectrum (expno) level references the default procpars.
"""
struct ParamDict{K, V} <: AbstractDict{K, V}
    acqupars::Dict{K, V}
    procpars::Dict{K, V}
    # Inner constructor
    function ParamDict(a::D, p::D) where D <: AbstractDict
        K = keytype(D)
        V = valtype(D)
        return new{K, V}(a, p)
    end
end
Base.getindex(d::ParamDict, key) = get(d.acqupars, key, get(d.procpars, key, throw(KeyError(key))))
Base.haskey(d::ParamDict, key) = haskey(d.acqupars, key) || haskey(d.procpars, key)

Base.iterate(d::ParamDict, state...) = iterate(merge(d.procpars, d.acqupars), state...)
Base.length(d::ParamDict) = length(merge(d.procpars, d.acqupars))
Base.keys(d::ParamDict) = keys(merge(d.procpars, d.acqupars))
Base.pairs(d::ParamDict) = pairs(merge(d.procpars, d.acqupars))

ParamDict(s::S) where S<:AbstractSpectrum = begin
    for k in keys(s.procs)
        s.procs[k].params = ParamDict(s.acqu, s.procs[k].params)
    end
end

#! TODO get this finished
""" 
    WrappedSpectrum{T, A} <: AbstractSpectrum
    - title
    - t : time domain data, fid/ser
    - s : frequency domain data (spectrum)
    - spectrum <: AbstractSpectrum 
    - f <: Function : processing functions mapping t -> s
    - ft <: AbstractFFT : f(t) -> φ(s)
------------------------------------------------------------------------
A wrapper which stores a raw spectrum file, and modified t and s domain
data. Initializes by reading in a fid and computing
    ψ(s) = fft(fftshift(f(t)))
with no other processing applied (note this includes SI zero filling).
For convenience, indexing operations are defined to access the frequency
domain array s by default, or if passed a string, to retrieve a stored 
parameter:
    s[i], s[i:j], s[(i,j,k...)], s[:]; i,j,k ∈ ℕ  -> s.s[... (idx)]
    s[a], s[a:b], s[(a,b,c...)];       a,b,c ∈ f64-> s.s[... (ppm)]
    s["par"]                                      -> s.spectrum["par"]
"""
mutable struct WrappedSpectrum{T <: AbstractFloat, A <: AbstractArray{T}} <: AbstractSpectrum
    title :: S where S <: AbstractString
    t :: A
    s :: A
    spectrum :: BrukerSpectrum
    plan :: p where p <: AbstractFFTs.Plan
end

# Convenience methods to access the spectrum
Base.real(s::WrappedSpectrum) = real(s.s)
Base.imag(s::WrappedSpectrum) = imag(s.s)
Base.abs(s::WrappedSpectrum) = abs(s.s)
Base.view(s::WrappedSpectrum, v) = @view s.s[v]
Base.size(s::WrappedSpectrum) = size(s.s)

# Indexing methods
Base.getindex(s::WrappedSpectrum, i::Int) = s.s[i]
Base.getindex(s::WrappedSpectrum, ::Colon) = s.s[:]
Base.getindex(s::WrappedSpectrum, a::AbstractArray) = s.s[a]
Base.getindex(s::WrappedSpectrum, rng::Tuple{Float64,Float64}) = s.s[ppmtoindex(s.spectrum,rng)]
Base.getindex(s::WrappedSpectrum, δ::Float64) = s.s[ppmtoindex(s.spectrum, δ)]
Base.getindex(s::WrappedSpectrum, param::String) = s.spectrum[param]

#! TODO
""" WrappedSpectrum(s::BrukerSpectrum)
------------------------------------------------------------------------
Outer constructor to wrap an expno, transforming the FID and storing FFT 
coeffs."""
function WrappedSpectrum(s::BrukerSpectrum)
    td = s["TD"]
    if length(td) == 1 
        t = s.fid
        dims = 1
    else
        # dynamically deduce the complex dims from FnMODE
        t = s.ser
        dims = (s["TD"], s["TD0"])
    end
end

Base.haskey(s :: BrukerSpectrum, k) = haskey(s.acqu, k) || haskey(s.procs[s.default_proc].params, k)

export BrukerSpectrum, PoptSpectrum, ProcessedSpectrum, AbstractSpectrum, WrappedSpectrum
