import Base.Threads: @threads, Atomic, atomic_add!
"""
    baseline_correct!(S, Δppm) -> S'
--------------------------------------------------------------------------
Normalises the intensity of the FT of each spectrum in an array to the 
average intensity of a region. 
    - S <: AbstractSpectrum
    - Δppm::Tuple{Float64,Float64} - should be a region of 𝑛𝑜𝑖𝑠𝑒.
"""
function baseline_correct!(s::S, ppm_range::Tuple) where S <: BrukerSpectrum
    ppm_range = float(ppm_range)
    idx1,idx2 = ppmtoindex.(ppm_range)
    s[:] .-= mean(s[idx1:idx2])
    s
end


function baseline_correct!(s::P, ppm_range::Tuple{Float64,Float64}) where P <: PoptSpectrum
    s[:] .-= mean(s[ppm_range])
    s
end


""" 
    intensity_normalise(s::Spectrum) -> (μ, σ) ∈ ℜ
--------------------------------------------------------------------------
Compute the appropriate normalising contstant for a given spectrum S. 
Works by measuring the the last 12.5% of the FID as noise, then accounting 
for the √n scaling w.r.t. n scans.
Returns μ describing the noise intensity of the spectrum and σ, the stdev 
of this value. 

Used for correcting for digitisation artefacts.
"""
function intensity_normalise(s :: S) where S <: AbstractSpectrum
    # get number of points as an integer    
    n = Int(s["TD"] // 8);
    μ = Atomic{Float64}(); μ[] = 0.;
    σ = Atomic{Float64}(); σ[] = 0.;
    fids = eachfibre(s.fid);
    idx = size(fids.arr, fids.dim) - n;
    #Threads.@threads 
    for fid in fids
        f = abs.(@view fid[idx:end]);
        μᵢ = mean(f);
        σᵢ = sqrt(sum((f .- μᵢ) .^ 2));
        atomic_add!(μ, μᵢ);
        atomic_add!(σ, σᵢ);
    end
    # noise scales as √n scans
    μ[] /= √n * fids.max
    # compute s.d. at the end, noting n/√n = √n
    σ[] /= (√n * fids.max) - 1
    return μ[], σ[]
end

export baseline_correct!, intensity_normalise
