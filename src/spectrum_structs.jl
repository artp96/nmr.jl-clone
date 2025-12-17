import Base.complex, Base.imag, Base.real
# Some utilities
const S = String
#Integration range in the frequency domain, values given in ppm."
const Intrng{T} = Tuple{T, T} where T <: AbstractFloat



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
    Outer Constructor definition.

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
struct BrukerSpectrum{A <: ℂᴺ, P <: ProcessedSpectrum}  <: AbstractSpectrum
    fid :: A
    acqu :: Dict{String, Any}
    procs :: Dict{Int, P}
    default_proc :: Int
    name :: String
    expno :: Int
        
end
    # Outer constructor to handle weakly typed dicts & split the "Real" FID.
function BrukerSpectrum(v :: V, a :: Dict{S, Any}, p :: Dict{I, P}, d :: I, n :: String, e :: I) where {
# Type-fu is intensifying #
        V <: ℝ¹, 
        P <: ProcessedSpectrum, 
        S <: AbstractString, 
        I <: Int}
        # process the FID, {ℝ, ℝ} -> ℂ
        v = split_fid(v)  # Ensure that the transformed FID is of the correct type
        if !(v isa AbstractVector{<: Complex})
            throw(ArgumentError("The split FID must be of type AbstractVector{$C}, but got $(typeof(f))"))
        end
        W = typeof(v)
        # Copy the acqupars into each procno
        s = BrukerSpectrum(v, a, p, d, n, e)
        ParamDict(s)
        return s
end

Base.getindex(s::BrukerSpectrum, i::Int) = s.procs[i]
Base.getindex(s::BrukerSpectrum, ::Colon) = s.procs[s.default_proc].re_ft
Base.getindex(s::BrukerSpectrum, a::AbstractArray) = s.procs[s.default_proc].re_ft[a]
Base.getindex(s::BrukerSpectrum, rng::Tuple{Float64,Float64}) = s[ppmtoindex(s,rng)]
Base.getindex(s::BrukerSpectrum, δ::Float64) = s[s.default_proc].re_ft[ppmtoindex(s, δ)]
Base.haskey(s :: BrukerSpectrum, k) = haskey(s.acqu, k) || haskey(s.procs[s.default_proc].params, k)
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

function BrukerSpectrum(fid :: V, acqu :: Dict{S, Any}, proc :: P) where 
    {
    T <: AbstractFloat, 
    V <: AbstractVector{T}, 
    S <: AbstractString, 
    P <: ProcessedSpectrum{T, V}
    }
    return BrukerSpectrum(fid, acqu, Dict{Int, ProcessedSpectrum}(1=>proc), 1, "", "")
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
struct PoptSpectrum{
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
        
        plan = plan_fft(ser, [1])
        
        
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
split_fid(fid :: V{ℝ}) where {V <: AbstractVector, ℝ<:Real} 
Split a fresh FID into it's real and complex parts.
N-D spectra are still acquired time-domain sequentially, so this should work for all N-D spectra.
Effectively applies the default zero-filling value of 2x in the fid dimension.
"""
function split_fid(fid :: V) where V <: ℝ¹
    fid = complex(fid)
    fid[1:2:end] .-= 1im * fid[2:2:end]
    fid[2:2:end] .-= fid[2:2:end]
    return fid
end
function split_fid(fid :: CV) where CV <: CuVector{<:f64}
    fid = complex(fid)
    fid .-= 1im * circshift(fid, -1)
    fid .-= circshift(fid, -1)
    return fid
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
import AbstractFFTs: AbstractFFTs.Plan
""" 
    WrappedSpectrum{T, A} <: AbstractSpectrum
    - title
    - t : time domain data, fid/ser
    - s : frequency domain data (spectrum)
    - src <: AbstractSpectrum 
    - fplan : FFTW plan for t -> s transform
    - iplan : FFTW plan for s -> t transform
-------------------------------------------------------------------------------------------
A wrapper which stores a raw spectrum file, and modified t and s domain
data. Initializes by reading in a fid from a BrukerSpectrum and applying
with no other processing applied (note this includes SI zero filling).
For convenience, indexing operations are defined to access the frequency
domain array s by default, or if passed a string, to retrieve a stored 
parameter. Modified parameters are stored in the params dict.
    s[i], s[i:j], s[(i,j,k...)], s[:]; i,j,k ∈ ℕ  -> s.src[... (idx)]
    s[a], s[a:b], s[(a,b,c...)];      a,b,c ∈ f64 -> s.src[... (ppm)]
    s["par"]                                      -> s.src["par"]
The FID is scaled by dividing by the noise level.
"""
mutable struct WrappedSpectrum{A<:ℂᴺ} <: AbstractSpectrum
    name :: S where S <: AbstractString
    ft :: A # time domain data
    fs :: A # freq domain data
    src :: BrukerSpectrum
    #fplan :: P where P <: AbstractFFTs.Plan
    #iplan :: P where P <: AbstractFFTs.Plan
    params :: D where D <: Dict{String, Any}
    expno :: Int64
    # Inner constructor method which re-gets dimstate
    function WrappedSpectrum(n::String,ft::A,fs::A,src::B,e::Int64) where {
    A<:ℂᴺ, B <: BrukerSpectrum
    }
        # get the structure of the FID - 2D, or pseudo-2D
        new{A}(n, ft, fs, src, Dict{String,Any}(), e)
    end
end

""" WrappedSpectrum(s::BrukerSpectrum)
-------------------------------------------------------------------------------------------
Outer constructor to wrap a Bruker expno, transforming the FID and storing FFT coeffs.
By default, stores on the GPU if available.
Scales the data by the noise level, assuming the last 12.5% of the FID is purely noise.
"""
function WrappedSpectrum(src::BrukerSpectrum)
    fs, ft = fft(src; give_processed_fid = true)
    return WrappedSpectrum(src.name, ft, fs, src, src.expno)
end

# Convenience methods to access the spectrum
Base.real(w::WrappedSpectrum) = real(w.fs)
Base.imag(w::WrappedSpectrum) = imag(w.fs)
Base.abs(w::WrappedSpectrum) = abs(w.fs)
Base.view(w::WrappedSpectrum, v) = @view w.fs[v]
Base.size(w::WrappedSpectrum) = size(w.fs)

# Indexing methods
Base.getindex(w::WrappedSpectrum, i::Int) = w.fs[i]
Base.getindex(w::WrappedSpectrum, ::Colon) = w.fs[:]
Base.getindex(w::WrappedSpectrum, a::AbstractArray) = w.fs[a]
Base.getindex(w::WrappedSpectrum, rng::Tuple{Float64,Float64}) = w.fs[ppmtoindex(w.src,rng)]
Base.getindex(w::WrappedSpectrum, δ::Float64) = w.fs[ppmtoindex(w.src, δ)]
Base.getindex(w::WrappedSpectrum, param::String) = begin 
    try 
        w.params[param]
    catch e;
        w.src[param];
    end
end

export BrukerSpectrum, PoptSpectrum, ProcessedSpectrum, AbstractSpectrum, WrappedSpectrum
