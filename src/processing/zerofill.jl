# Contains a function to apply zero filling to an FID and several convenience methods to 
# apply the default zero filling implied by SI (Bruker).

export zero_fill

""" zero_fill(v :: Vector, k :: Int) -> v[1...end 0 0 0 ...0] 
-------------------------------------------------------------------------------------------
Bottom level function to append k zeroes to the FID stored in vector v, in place.
"""
function zero_fill(v :: V, k :: Int64) where V <: AbstractVector
    k < length(v) && @warn "k < TD, truncating FID by $(-k+s["TD"])!"
    sizehint!(v, length(v) + k)
    append!(v, zeros(eltype(v), k))
    return v
end

""" zero_fill( s <: BrukerSpectrum) -> [s.fid,0 0 0...]
--------------------------------------------------------------------------------------------
By default, applies the ZF implied by the spectral `params`.
"""
zero_fill(s :: BrukerSpectrum) = zero_fill(s.fid, floor(Int, s["SI"] - s["TD"]/2))
