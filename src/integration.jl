export intensity, integrate
intrng(s::BrukerSpectrum) = s[s.default_proc].intrng
intrng(p::PoptSpectrum) = p.proc.intrng
intrng(w::FracSpectrum) = intrng(w.src)

function intrng_indices(s::BrukerSpectrum)
    rng = intrng(s)
    [ppmtoindex(s,i[1]):ppmtoindex(s,i[2]) for i in rng]
end

function intrng_data(s::BrukerSpectrum)
    [@view s[r] for r in intrng_indices(s)]
end

function intrng_shifts(s::BrukerSpectrum)
    rng = intrng(s)
    [range(i[1]; stop=i[2], length=ppmtoindex(s,i[2])-ppmtoindex(s,i[1])+1) for i in rng]
end

function find_rng(p::ProcessedSpectrum, δ::Float64)
    for (i,(hi,lo)) in enumerate(p.intrng)
        if lo < δ < hi
            return i
        end
    end
end

find_rng(s::BrukerSpectrum, δ::Float64) = find_rng(s[s.default_proc], δ)

function remove_rng(sp::Union{BrukerSpectrum,ProcessedSpectrum}, δ::Float64)
    t = deepcopy(sp)
    remove_rng!(t, δ)
    t
end

remove_rng!(s::BrukerSpectrum, δ::Float64) = remove_rng!(s[s.default_proc], δ)
remove_rng!(p::ProcessedSpectrum, δ::Float64) = remove_rng!(p, find_rng(p, δ))
remove_rng!(p::ProcessedSpectrum, n::Int) = deleteat!(p.intrng, n)
remove_rng!(p::ProcessedSpectrum, x::Nothing) = nothing

integrate(v::Vector, r::UnitRange) = sum(v[r])
integrate(p::ProcessedSpectrum, r::UnitRange) = integrate(p[:], r)
integrate(s::BrukerSpectrum, r::UnitRange) = integrate(s[s.default_proc], r)

function integrate(s::BrukerSpectrum, ppm_range::Tuple{Float64,Float64})
    r1 = ppmtoindex(s, ppm_range[1])
    r2 = ppmtoindex(s, ppm_range[2])
    integrate(s, r1:r2)
end

function integrate(s::BrukerSpectrum)
    [integrate(s, r) for r in intrng(s)]
end

function integrate(s::BrukerSpectrum, ref_rng::Int)
    rngs = intrng(s)
    ref_int = integrate(s, rngs[ref_rng])
    integrate(s)./ref_int
end

integrate(s::BrukerSpectrum, δ::Float64) = integrate(s, find_rng(s, δ))

#integrate(w::FracSpectrum, rng::UnitRange, kind = real) = 

function integrate(w::FracSpectrum, ppm_range::Tuple{Float64,Float64}, kind = real)
    r1 = ppmtoindex(w, ppm_range[1])
    r2 = ppmtoindex(w, ppm_range[2])
    return sum(kind(w)[r1:r2])
end
###########################################################################


intensity(v::Vector, r::UnitRange) = maximum(v[r])
intensity(p::ProcessedSpectrum, r::UnitRange) = intensity(p[:], r)
intensity(s::BrukerSpectrum, r::UnitRange) = intensity(s[s.default_proc], r)

function intensity(s::BrukerSpectrum, ppm_range::Tuple{Float64,Float64})
    r1 = ppmtoindex(s, ppm_range[1])
    r2 = ppmtoindex(s, ppm_range[2])
    intensity(s, r1:r2)
end

function intensity(s::BrukerSpectrum)
    [intensity(s, r) for r in intrng(s)]
end

function intensity(s::BrukerSpectrum, ref_rng::Int)
    rngs = intrng(s)
    ref_int = intensity(s, rngs[ref_rng])
    intensity(s)./ref_int
end

intensity(s::BrukerSpectrum, δ::Float64) = intensity(s, find_rng(s, δ))

function intensity(w::FracSpectrum, ppm_range::Tuple{Float64,Float64}, kind = real)
    r1 = ppmtoindex(w, ppm_range[1])
    r2 = ppmtoindex(w, ppm_range[2])
    return maximum(kind(w)[r1:r2])
end
