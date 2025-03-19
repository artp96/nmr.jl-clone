

# Some utilities
"Integration range in the frequency domain, values given in ppm."
const Intrng{T} = Tuple{T, T} where T <: Union{AbstractFloat, Missing}
abstract type AbstractSpectrum end

"""
Processed (frequency domain) NMR spectrum. Parameters follow the Bruker conventions:
- **re_ft**: real part of Fourier transform
- **im_ft**: imaginary part of Fourier transform
- **params**: dictionary of processing parameters, mostly from `proc` file for Bruker data
- **intrng**: list of integration ranges
- **title**: Spectrum title
"""

mutable struct ProcessedSpectrum{T <: AbstractFloat, V <: AbstractVector{T}, S <: AbstractString} <: AbstractSpectrum
    re_ft :: V
    im_ft :: V
    params :: Dict{S, Any}
    intrng :: Intrng
    procno :: Int
    title :: String
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
mutable struct Spectrum{T <: AbstractFloat, V <: AbstractVector{T}, P <: ProcessedSpectrum{T, V}}  <: AbstractSpectrum
    fid :: V
    acqu :: Dict{String, Any}
    procs :: Dict{Int, P}
    default_proc :: Int
    name :: String
    expno :: Int
        
        """Inner constructor to handle weakly typed dicts """
    function Spectrum(f :: V, a :: Dict{S, Any}, p :: Dict{I, P}, d :: I, n :: String, e :: I) where {
        # Type-fu is intensifying #
        T <: AbstractFloat, V <: AbstractVector{T}, P <: ProcessedSpectrum{T, V}, S <: AbstractString, I <: Int}
            
        # for (k, v) in (collect(keys(a)), values(a))
        #     a_new = Dict{String, Any}(string(k) => v)
        # end
        # for (k, v) in (collect(keys(p)), values(p))
        #     p_new = Dict{Int, ProcessedSpectrum{T, V}}(k => v)
        # end
        return new{T, V, P}(f, a, p, d, n, e)
    end
end

mutable struct PoptSpectrum{T <: AbstractFloat, V <: AbstractVector{T}, A <: AbstractArray{T}} <: AbstractSpectrum
    fid :: V
    acqu :: Dict{String, Any}
    procs :: Dict{Int,ProcessedSpectrum{T, V}}
    default_proc :: Int
    name :: String
    expno :: Int
    popt_serfile :: A

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

Spectrum(fid :: V, acqu :: Dict{S, Any}, proc :: P) where {
    T <: AbstractFloat,
    V <: AbstractVector{T},
    S <: AbstractString,
    P <: ProcessedSpectrum{T, V}
    } = Spectrum(fid, acqu, Dict{Int, ProcessedSpectrum{T, V}}(1=>proc), 1, "", "")

export Spectrum
