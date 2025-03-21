# Made with GPT-4o
"""
    FiberIterator(arr::AbstractArray, dim::Int)

Creates an iterator that returns each 1D fiber of `arr` along the specified dimension `dim`.
"""
struct FiberIterator{T, N}
    arr::AbstractArray{T, N}
    dim::Int
end

# Make FiberIterator an iterator
function Base.iterate(fib::FiberIterator, state=(1,))
    idx = state[1]
    # Size of the array along the desired dimension
    size_dim = size(fib.arr, fib.dim)

    # If the index is out of range, stop iteration
    if idx > size_dim
        return nothing
    end

    # Create a view for the current fiber along the chosen dimension
    fiber = view(fib.arr, Base.OneTo(size(fib.arr, 1))..., idx)
    
    # Move to the next index
    return fiber, (idx + 1,)
end

function Base.IteratorSize(::Type{<:FiberIterator})
    Base.HasLength()
end

function Base.length(fib::FiberIterator)
    size(fib.arr, fib.dim)
end

# slighlty simpler interface with default behaviour matching 
# the structure of popt serfiles
eachfibre(a; dim = 1) = FiberIterator(a, dim)
export eachfibre
