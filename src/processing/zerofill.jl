# Contains a function to apply zero filling to an FID and several convenience methods to 
# apply the default zero filling implied by SI (Bruker).
export zero_fill
import Base: ndims 
import SparseArrays: sparse!
ndims(idx::CartesianIndex) = 1

""" zero_fill(v :: Vector, k :: Int) -> v[1...end 0 0 0 ...0] 
-------------------------------------------------------------------------------------------
Bottom level function to append zeroes to the FID stored in vector v, until it has length k
"""
function zero_fill(v :: V, k :: Int64) where V <: AbstractVector
    if k < length(v) 
        #@warn "k < TD, truncating FID by $(length(v)-k)!"
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

""" zero_fill(A, k) A∈ℂ¹ʲ -> [A₁ (0 0 ...)ₖ A₂ (0 0 ...)ₖ ...] ∈ 𝕊¹ᵏʲ 
-------------------------------------------------------------------------------------------
Bottom level function to zero fill an array in a given dimension by promoting its sparsity.
Called on pre-zero-filled sparse arrays, alters the sparsity to the new k (not cumulative).
"""
zero_fill(v::𝕍, n::n64) where 𝕍<:ℂ¹ = n < 2 ? v : sparsevec(n*eachindex(v) .- 1, v)
zero_fill(A::𝔸, facs::N) where {𝔸<:ℂ², N<:ℕ¹} = begin
    length(facs) != ndims(A) && throw(ArgumentError("Calls to zero_fill arrays should include a scaling factor ≥1 in each dimension."))
    # get the current idxs
    dense_idxs = CartesianIndices(A)
    sparse_idxs = similar(dense_idxs)
    # apply the scaling factor in each dimension
    [sparse_idxs[i] = CartesianIndex((dense_idxs[i].I .* facs)...) for i ∈ eachindex(dense_idxs)]
    # reshape the sparsified indexes into vectors in each dimension
    sparse_idxs = Tuple([[getindex.(sparse_idxs, d)...] .- (facs[d] - 1) for d in 1:ndims(sparse_idxs)])

    # return the sparsified array with fac structural zeros in each dimension between entries, and after them.
    return sparse!(sparse_idxs..., [A[i] for i in eachindex(A)], facs .* size(A)...)

end

# calling zerofill on a sparse array should override the previous zerofilling pattern
zero_fill(sparse::M, facs::N) where {M<:𝕊², N<:ℕ¹} = begin
    length(facs) != ndims(A) && throw(ArgumentError("Calls to zero_fill arrays should include a scaling factor ≥1 in each dimension."))
    I, J, V = findnz(sparse)
    # deduce the old zero filling factors, and divide the new ones
    dI = I[findfirst(i -> i ≠ I[1])] - I[1]
    dJ = J[findfirst(j -> j ≠ J[1])] - J[1]
    facs ./= (dI, dJ)
    [sparse_idxs[v] = CartesianIndex(((I[v], J[v]) .* facs)...) for v ∈ eachindex(V)]
    sparse_idxs = Tuple([[getindex.(sparse_idxs, d)...] .- (facs[d] - 1) for d in 1:ndims(sparse_idxs)])
    return sparse!(sparse_idxs..., V, facs .* size(A)...)
end
