# Contains a function to apply zero filling to an FID and several convenience methods to 
# apply the default zero filling implied by SI (Bruker).
export zero_fill
import Base: ndims 
import SparseArrays: sparse!
ndims(idx::CartesianIndex) = 1

#=
""" zero_fill( s <: BrukerSpectrum) -> [s.fid,0 0 0...]
--------------------------------------------------------------------------------------------
By default, applies the ZF implied by the spectral `params`. If SI < TD, by default returns 
the full FID.
"""
zero_fill(s :: BrukerSpectrum) = zero_fill(s.fid, s["SI"]÷length(fid))

""" zero_fill(A, k) A∈ℂ¹ʲ -> [A₁ (0 0 ...)ₖ A₂ (0 0 ...)ₖ ...] ∈ 𝕊¹ᵏʲ 
-------------------------------------------------------------------------------------------
Bottom level function to zero fill an array in a given dimension by promoting its sparsity.
Called on pre-zero-filled sparse arrays, alters the sparsity to the new k (not cumulative).
"""
zero_fill(v::𝕍, n::n64) where 𝕍<:ℂ¹ = n > 1 ? sparsevec((n*eachindex(v)) .- (n-1), v, n*length(v)) : @error "zero filling is multiplicative, n = $n < 0!"

function zero_fill(m::M, facs::N) where {M<:Union{ℝ², ℂ²}, N<:ℕ¹}
    length(facs) != ndims(m) && throw(ArgumentError("Calls to zero_fill arrays should include a scaling factor ≥1 in each dimension."))
    # get the current idxs
    dense_idxs = CartesianIndices(m)
    sparse_idxs = similar(dense_idxs)
    # apply the scaling factor in each dimension
    [sparse_idxs[i] = CartesianIndex((dense_idxs[i].I .* facs)...) for i ∈ eachindex(dense_idxs)]
    # reshape the sparsified indexes into vectors in each dimension
    sparse_idxs = Tuple([[getindex.(sparse_idxs, d)...] .- (facs[d] - 1) for d in 1:ndims(sparse_idxs)])

    # return the sparsified array with fac structural zeros in each dimension between entries, and after them.
    return sparse!(sparse_idxs..., [m[i] for i in eachindex(m)], facs .* size(m)...)
end

# calling zerofill on a sparse array should override the previous zerofilling pattern
zero_fill(sparse::M, facs::N) where {M<:𝕊², N<:ℕ¹} = begin
    length(facs) != ndims(A) && throw(ArgumentError("Calls to zero_fill arrays should include a scaling factor ≥1 in each dimension."))
    I, J, V = findnz(sparse)
    # deduce the old zero filling factors, and divide the new ones
    dI = I[findfirst(i -> i ≠ I[1])] - I[1]
    dJ = J[findfirst(j -> j ≠ J[1])] - J[1]
    facs ./= (dI, dJ)
    sparse_idxs = similar(CartesianIndex.(I,J))
    [sparse_idxs[v] = CartesianIndex(((I[v], J[v]) .* facs)...) for v ∈ eachindex(V)]
    sparse_idxs = Tuple( [[getindex.(sparse_idxs, d)...] .- (facs[d] - 1) for d in 1:ndims(sparse_idxs)] )
    return sparse!(sparse_idxs..., V, facs .* size(sparse)...)
end

# zero filling CuVector or CuMatrixes using sparsity is rarely worth it as the followup calculation is generally O(C). 
# index manipulation is much more easily done on the CPU in this case, especially as this method avoids allocations 
# of zeros on the GPU where memory is more scarce.
zero_fill(m::CuA, args::n64) where {CuA<:CuArray} = zero_fill(Array(m), args) |> cu
zero_fill(m::CuA, args::ℕ¹) where {CuA<:CuArray} = zero_fill(Array(m), args) |> cu
#!TODO this method needs another look, it could get expensive, and is very much the quick-and-dirty version
zero_fill(m::CuS, args::n64) where CuS<:CUSPARSE.AbstractCuSparseArray = zero_fill(sparse(Array(m)), args)
zero_fill(m::CuS, args::ℕ¹) where CuS<:CUSPARSE.AbstractCuSparseArray = zero_fill(sparse(Array(m)), args)


# needed for specialised method extension
import SparseArrays: sparsevec
"""
    Internal specialised method to correctly sparsify CuVectors
"""
sparsevec(idxs::ℕ¹, v::T, n::n64) where {Tv<:CuArray, Tt<:Tuple, T<:SubArray{c64, 1, Tv, Tt, true}} = CUDA.CUSPARSE.CuSparseVector(gpu!(idxs), CUDA.zeros(c64,length(v)).+v, n)

=#
#! TODO:
# This is not the way to zero fill spectra. In principle, zero_filling should be *free* if it has repeating structure.
# A prototype better solution is below.
#
"""
    AbstractArray wrapper which includes zero_filling factors in each dimension.
    Where only a single factor is passed, it is applied to every dimension.
"""
struct FID{N,A<:AbstractArray{c64,N}} <: AbstractArray{c64,N}
   data :: A
   facs :: NTuple{N,Int}
   function FID(data::A, facs::T) where {N, A<:AbstractArray{c64,N},T<:NTuple}
       @assert length(facs) == N "One factor per dimension required"
       @assert all(>(0), facs) "All factors must be positive" 
       new{N, A}(data, facs)
   end
"""
    Gpu constructor explicitly stores zeros.
"""
    function FID(data::X, facs::T) where {N, X<:CuArray{c64,N},T<:NTuple}
        outdim = size(data) .* facs
        expanded = CUDA.zeros(c64, outdim)
        strides = ntuple(d -> stride(expanded, d), N)
        nt = 1024
        @cuda threads=nt blocks=cld(length(data),nt) _gpu!_FID_scatter_kernel!(expanded, data, facs, strides)
        return new{N, X}(expanded, facs)
    end
end
function _gpu!_FID_scatter_kernel!(out, in, facs, strides)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    i > length(in) && return
    I = CartesianIndices(in)[i]
    # Initialise the target index in the output array.
    out_idx = 1
    @inbounds for d in 1:length(facs)
        out_idx += ((I[d] - 1) * facs[d]) * strides[d]
    end

    @inbounds out[out_idx] = in[i]
    return nothing
end
"""
    Move an FID to the GPU, storing only the nonzero values.
"""
gpu!(f::FID{N, A}) where {N, A<:Array} = FID(gpu!(f.data), f.facs)
"""
    Move an FID to the CPU, storing only the nonzero values.
"""
cpu!(f::FID{N, X}) where {N, X<:CuArray} = begin
    idxs = CartesianIndices(f)
    data = @view f.data[idxs]
    data = collect(data)
    return FID(data, f.facs)
end
Base.IndexStyle(::Type{<:FID}) = IndexCartesian()
Base.size(f::FID{N,A}) where N where A<:Array = size(f.data) .* f.facs
Base.size(f::FID{N,X}) where N where X<:CuArray = size(f.data)
Base.show(io::IO, fid::FID) = begin
    show(io, "FID with zero-filling factor $(fid.facs) in each dimension. Data:\n")
    show(io, fid.data) 
end
Base.display(fid::FID) = begin
    display("FID with zero-filling factor $(fid.facs) in each dimension. Data:")
    display(fid.data) 
end
Base.length(f::FID) = prod(size(f))
Base.CartesianIndices(f::FID{N}) where N = CartesianIndices(ntuple(d -> 1:f.facs[d]:size(f)[d], N))
nnz(f::FID) = length(f.data)
# Helper: check if index lies on the grid (1-based)
@inline function _on_grid(i::Int, fac::Int)
    ((i - 1) % fac) == 0
end

@inline function _map_index(i::Int, fac::Int)
    ((i - 1) ÷ fac) + 1
end

function Base.getindex(f::FID{N}, I::Vararg{Int,N}) where {N}
    @inbounds begin
        for d in 1:N
            if !_on_grid(I[d], f.facs[d])
                return zero(c64)
            end
        end
        J = ntuple(d -> _map_index(I[d], f.facs[d]), N)
        return f.data[J...]
    end
end

function Base.setindex!(f::FID{N}, x::c64, I::Vararg{Int,N}) where {N}
    @inbounds begin
        for d in 1:N
            if !_on_grid(I[d], f.facs[d])
                return f
            end
        end
        J = ntuple(d -> _map_index(I[d], f.facs[d]), N)
        f.data[J...] = x
        return f
    end
end

"""
    Lazy cpu-side constructor.
"""
FID(a::A) where A <: AbstractArray = FID(a, ntuple(_ -> 1, ndims(a)))
zero_fill(a::A, t) where A<:AbstractArray = FID(a, t)
"""
    Filling an n-D array with one factor fills the entire array.
"""
zero_fill(data::A, facs::Int) where {N,A<:AbstractArray{c64}} = FID(data, ntuple(_ -> facs, ndims(A)))
# zero_fill should overwrite the previous filling factors.
zero_fill(f::FID{N, A}, facs::T) where {N,A<:Array{c64, N}, T<:NTuple{N}} = FID(f.data, facs)
# Reindexing is done most efficiently on the CPU.
zero_fill(f::FID{N, X}, facs::T) where {N,X<:CuArray{c64, N}, T<:NTuple{N}} = begin
    idxs = CartesianIndices(f)
    data = @view f.data[idxs]
    data = collect(data)
    return FID(gpu!(data), facs)
end


# """
#     Convenience constructor: scalar factor applies to all dimensions
# """
# FID(data::A, fac::Int) where {A<:AbstractArray{c64}} =
#     FID(data, ntuple(_ -> fac, ndims(data)))

x_1 = rand(c64, 10)
x_2 = zero_fill(x_1, 2)

X_1 = rand(c64, 10, 10)
X_2 = zero_fill(X_1, 2)


g_1 = CUDA.rand(c64, 10)
g_2 = zero_fill(g_1, 2)

G_1 = CUDA.rand(c64, 10, 10)
G_2 = zero_fill(G_1, 2)
