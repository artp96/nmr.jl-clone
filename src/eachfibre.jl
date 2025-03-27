# Made with GPT-4o's help, I don't get iterators really.
"""
    fibreIterator(arr::AbstractArray, dim::Int)

Creates an iterator that returns each 1D fiber of `arr` along the specified dimension `dim`.
Like eachcol and eachrow, this is effectively a wrapper for Slices.
"""
struct fibreIterator{T, N}
    arr :: AbstractArray{T, N}
    dim :: Int
    max :: Int
    fibres :: Slices
    function fibreIterator(arr :: AbstractArray{T, N}, dim::Int) where {T, N} 
        any(size(arr) .== 0) && Throw(ArgumentError("At least one dimension of $arr is empty! \nsize = $(size(arr))."))
        max = prod(size(arr)) ÷ size(arr, dim)
        slice_dims = filter( d -> d != dim, ntuple(identity, ndims(arr)) )
        fibres = eachslice(arr; dims = slice_dims)
        new{T, N}(arr, dim, max, fibres)
    end
end

# Make fibreIterator an iterator
function Base.iterate(fib::fibreIterator, state=(1,))
    idx = state[1]

    # If the index is out of range, stop iteration
    if idx > fib.max
        return nothing
    end

    # views do not work with CUDA, see cuda.subview()
    # Move to the next index
    return fib.fibres[idx], (idx + 1,)
end

function Base.IteratorSize(::Type{<:fibreIterator{T, N}}) where {T, N <: Integer}
    M = N - 1
    Base.HasShape{M}()
end

function Base.getindex(fib::fibreIterator, i::Int)
    fib.fibres[i]
end


Base.size(fib::fibreIterator) = size(fib.fibres)

function Base.getindex(fib::fibreIterator, inds...)

    fib.fibres[i]
end
function Base.length(fib::fibreIterator)
    return fib.max
end

Base.firstindex(fib::fibreIterator) = !isempty(fib) ? 1 : Throw(ArgumentError("fibreIterator is empty!"))
Base.lastindex(fib::fibreIterator) = !isempty(fib) ? fib.max : Throw(ArgumentError("fibreIterator is empty!"))

# slighlty simpler interface with default behaviour matching 
# the structure of popt serfiles
eachfibre(a :: A; dim = 1) where {A <: AbstractArray} = fibreIterator(a, dim)

export eachfibre, fibreIterator
