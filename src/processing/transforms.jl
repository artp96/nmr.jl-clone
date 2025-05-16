#=
This file contains some specialised versions of Fourier transforms to handle spectrum data 
structures.
=#

using AbstractFFTs: Plan
import FFTW: fft, ifft, plan_fft,  plan_ifft

todo = !ismissing("yes!")

""" fft(s <: AbstractSpectrum; f(t) = s.fid) -> ψ ∈ ℂᴺ.
-------------------------------------------------------------------------------------------
Forward Fourier transform of an spectrum which calls a specialsed method depending on the 
spectrum in question. Returns a transformed array, by default phase corrected using the 
Bruker parameters PHC0 and PHC1 respectively.
Optionally, pass a specific array as fid to test the processing associated with 's'. 
"""
fft(s::S) where S<:AbstractSpectrum = fft(s, s["PHC0"], s["PHC1"])
fft(w::W; ft::A=w.ft) where {W<:WrappedSpectrum,A<:AbstractArray} = fft(w, w["PHC0"], w["PHC1"]; ft=ft)

function fft(w::WrappedSpectrum, ϕ₀::T, ϕ₁::T; ft=w.ft) where 
    {
    T<:AbstractFloat,
    } 
    k = get_DSF_offset(w.src)
    # If there is a residual first-order phase offset from the DSP shift, add this to the 
    # given ϕ₁.
    if !iszero(k % 1)
        ϕ₁ += (k % 1) * 360 * size(w.ft, 1) / ( w.src["SW"] * w.src["SFO1"] )
        k = floor(Int, k)
    end
    # transform the spectrum efficiently using the fft plan
    fs = fftshift(fft(ft), 1)
    # Apply the phase correction in f1, only.
    for fₖs in eachfibre(fs)
        ϕ_correct!(fₖs, ϕ₀, ϕ₁)
    end
    # return the transformed FID.
    return fs
end

# fft the FID of a Bruker Spectrum, returning the spectrum and plan used. Applies
# default zero filling based on SI and phase correction based on PHC0/1.
function fft(s::BrukerSpectrum, ϕ₀::T=s["PHC0"], ϕ₁::T=s["PHC1"]; give_processed_fid = false) where
    {
        T<:AbstractFloat,
    } 
    ft = s.fid
    if isempty(ft) & give_processed_fid 
        @warn "Spectrum $(s.expno) - No fid to transform. Inverting the default processed data."
        return complex(s.procs[s.default_proc]), ifft(s.procs[s.default_proc])
    elseif isempty(ft) !give_processed_fid
        @warn "Spectrum $(s.expno) - No FID to transform. Returning the default processed data."
        return complex(s.procs[s.default_proc])
    end
        
    dims = s["FnMODE"] ≥ 2 ? [1, 2] : [1]
    k = get_DSF_offset(s)    
    chunk = round(Int, s["TD"] / 8)
    # If there is a residual first-order phase offset from the DSP shift, add this to the 
    # given ϕ₁.
    if !iszero(k % 1)
        ϕ₁ += (k % 1) * 360 * s["SI"] / ( 2s["SW"] * s["SFO1"] )
        k = floor(Int, k)
    end
    for fₖt in eachfibre(ft)
        # centre the noise at zero, ie. correct for a baseline offset, iff. ∃
        fₖt .-= mean(fₖt[chunk:end])

        # scale by the noise-per-scan
        # this should correct for some DSP artefacts so that signal correctly 
        # scales as √N scans.
        noise = mean(abs, fₖt[chunk:end]) / √s["NS"]
        fₖt ./= noise
    end
    # apply the default zero filling
    ft = zero_fill(s) |> gpu
    # fourier transform the FID
    fs = circshift(ft, k)
    fs = fft(fs, dims)
    fs = fftshift(fs, 1)
    # Apply the phase correction in f1, only.
    for fₖs in eachfibre(fs)
        ϕ_correct!(fₖs, ϕ₀, ϕ₁)
    end
    # return the transformed ft.
    if give_processed_fid
        return fs, ft
    else
        return fs
    end
end

# Fourier transform a fid and apply an analytical first order phase offset 
# to the t₁ dimension in 1 or n-D spectra, based on DSP. 
# ϕ₀ in °, ϕ₁ in ° Hz⁻¹.
#=
function fft(fid :: A, ϕ₀::T = 0.0, ϕ₁::T = 0.0, plan :: P = plan_fft(fid)) where {
    T <: AbstractFloat,
    A <: AbstractArray,
    P <: AbstractFFTs.Plan,
}
end
=# # general formulation may be redundant?  

# Bottom-level iffts for containers which hold frequency-domain spectra.
""" ifft(p <: ProcessedSpectrum) -> f(t)
-------------------------------------------------------------------------------------------
Generate an FID from the real and imaginary parts of a processed Bruker spectrum.
""" 
ifft(p::P) where P<:ProcessedSpectrum = ifft(ifftshift(complex(p)))


""" ifft(W <: WrappedSpectrum) -> f(t)
-------------------------------------------------------------------------------------------
Generate an FID from the real and imaginary parts of an NMR.jl-wrapped spectrum.
"""
ifft(w::WrappedSpectrum) = ifft(ifftshift(w.s))

#! TODO fix this / have stateful transforms
"""dimstate(s::AbstractSpectrum) -> [Int]
-------------------------------------------------------------------------------------------
Determine the dimensionality of the fourier transform associated with the spectrum S from 
FnMODE. 
"""

