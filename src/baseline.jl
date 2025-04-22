"""
    baseline_correct!(S, Δppm) -> S'
--------------------------------------------------------------------------
Normalises the intensity of the FT of each spectrum in an array to the 
average intensity of a region. 
    - S <: AbstractSpectrum
    - Δppm::Tuple{Float64,Float64} - should be a region of 𝑛𝑜𝑖𝑠𝑒.
"""
function baseline_correct!(s::S, ppm_range::Tuple{Float64,Float64}) where S <: Spectrum
    s[:] .-= mean(s[ppm_range])
    s
end


function baseline_correct!(s::P, ppm_range::Tuple{Float64,Float64}) where P <: PoptSpectrum
    s[:] .-= mean(s[ppm_range])
    s
end
