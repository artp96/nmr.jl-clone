# Contains a function to apply zero filling to an FID and several convenience methods to 
# apply the default zero filling implied by SI (Bruker).

export zero_fill

""" zero_fill(v :: Vector, k :: Int) -> v[1...end 0 0 0 ...0] 
-------------------------------------------------------------------------------------------
Bottom level function to append zeroes to the FID stored in vector v, until it has length k
"""
function zero_fill(v :: V, k :: Int64) where V <: AbstractVector
    if k < length(v) 
        @warn "k < TD, truncating FID by $(length(v)-k)!"
        return w = v[1:k]
    elseif k > length(v)
        w = v
        sizehint!(w, k)
        append!(w, zeros(eltype(v), k-length(v)))
    return w
    else
        return v
    end
end

""" zero_fill( s <: BrukerSpectrum) -> [s.fid,0 0 0...]
--------------------------------------------------------------------------------------------
By default, applies the ZF implied by the spectral `params`. If SI < TD, by default returns 
the full FID.
"""
zero_fill(s :: BrukerSpectrum) = zero_fill(s.fid, s["SI"])
