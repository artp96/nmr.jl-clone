# This script contains functions to map gradient nonuniformity across a 
# pseudo-2D spectrum.
#
# TODO:
# 1 - gradient nonuniformity mapping function 
# 2 - NUG Correction function 
# 3 - B1-NUG mapping
    
map_gradient_nonuniformity(s::BrukerSpectrum; kwargs...) = begin
    w = FracSpectrum(s)
    return map_gradient_nonuniformity(w; kwargs...)
end



function map_gradient_nonuniformity(w::FracSpectrum)
    # unpack (shallowly) the variables we need
    arr = w.s;
    gpz = w

end


# product of gaussian signal profile with bandwidth w centred at spoffs s₀, 
# a linear gradient signal profile, and a d²g/dz² term
#profile(w, s, s₀ = 0.,) = exp(-( (s-s0)/2.355w )^2 ) * ( (s - s₀) + (s^3 - s₀) )
""" freq_profile(s, g)
    - s <: BrukerSpectrum
    - g - applied gradient, default "GPZ9" * GMAX
--------------------------------------------------------------------
"""
function freq_profile(s, g)
    s_ζ = real(s)    
    
end
