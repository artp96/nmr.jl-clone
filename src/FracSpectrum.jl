""" 
    FracSpectrum{A, N, B} <: AbstractSpectrum
    - title : string
    - arr <: AbstractArray{ℂ, N} : spectrum data.
    - a <: Tuple{f64, N}: transform angle, in π / 2 radians, of each dimension.
    - src :: B <: AbstractSpectrum
    - params :: D where D <: Dict{String, Any}
-------------------------------------------------------------------------------------------
Called to wrap an imported spectrum. Tracks the time-frequency state of the spectrum in a. 
Each dimension can be transformed independently.
The stored time-domain data is always at least preprocessed, ready for processing and 
transformation. For example, 
    - It has the correct dimensions.
    - A Bruker FID has been made complex, and had it's GRPDLY offset corrected.
    - A Bruker Procno spectrum has been combined and made complex, unless it is real only.
Note that α is given in units of radians. According to convention, this variable could be 
called "a" and be given in units of π/2 rad. such that "a = 1" corresponds to the Fourier
domain and a = 0 to the time domain.
See also the convenience functions α(a) and get_u(a).
"""
mutable struct FracSpectrum{N, B} <: AbstractSpectrum where {N, B<:AbstractSpectrum}
    name :: S where S <: AbstractString
    s :: A where A <: ℂᴺ{N} # Spectrum data.
    src :: B
    α :: NTuple{N, f64}
    params :: D where D <: Dict{String, Any}
    # Inner constructor method which re-gets dimstate
function FracSpectrum(n::String, s::A, src::B, α, p=Dict{String,Any}) where {N,A<:ℂᴺ{N},B<:AbstractSpectrum}
        @assert length(α) == ndims(s) "Must have an α value for each dimension of the spectrum."
        # get the structure of the FID - 2D, or pseudo-2D
        new{N,B}(n, gpu!(s), gpu!(src), ntuple(i->float(α[i]),N), p)
    end
end

""" FracSpectrum(s::BrukerSpectrum)
-------------------------------------------------------------------------------------------
Outer constructor to wrap a Bruker expno, transforming the FID and storing FFT coeffs.
By default, stores on the GPU if available.
Scales the data by the noise level, assuming the last 12.5% of the FID is purely noise.
"""
function FracSpectrum(src::BrukerSpectrum)
    s = preprocess(src)
    α = zeros(f64, ndims(s))
    return FracSpectrum(src.name, s, src, α)
end

# Convenience methods to access the spectrum
Base.real(w::FracSpectrum) = real(w.s)
Base.imag(w::FracSpectrum) = imag(w.s)
Base.abs(w::FracSpectrum) = abs(w.s)
Base.view(w::FracSpectrum, v) = @view w.s[v]
Base.size(w::FracSpectrum) = size(w.s)

# Indexing methods
Base.getindex(w::FracSpectrum, i::Int) = w.s[i]
Base.getindex(w::FracSpectrum, ::Colon) = w.s[:]
Base.getindex(w::FracSpectrum, a::AbstractArray) = w.s[a]
Base.getindex(w::FracSpectrum, rng::Tuple{Float64,Float64}) = w.s[ppmtoindex(w.src,rng)]
Base.getindex(w::FracSpectrum, δ::Float64) = w.s[ppmtoindex(w.src, δ)]
Base.getindex(w::FracSpectrum, param::String) = begin 
    try 
        w.params[param]
    catch e;
        # iterates into procno params if appropriate.
        w.src[param];
    end
end


