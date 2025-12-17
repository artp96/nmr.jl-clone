# Dynamically deduce the presence of a functional CUDA gpu and set the 
# gpu interface up. In-place definitions ensure the objects only exist in one place.
dense!(x::S) where S<:𝕊ᴺ = return x = Array(x)
if has_cuda_gpu() & has_cuda()
    @info "CUDA found. GPU mode, use gpu!(data) and cpu!(data) to move data."
    """
        Move an array to the gpu!. Alias for CuArray, preserves sparsity of sparse complex arrays.
    """
    gpu!(x::A) where A<:AbstractArray = return x = CuArray(x)
    """
        Move an array to the cpu!. Alias for Array, preserves sparsity of sparse complex arrays.
    """
    cpu!(x::A) where A<:AbstractArray = return x = Array(x)
    # Automatically preserve sparsity.
    gpu!(x::S) where S<:𝕊¹ = return x = CUDA.CUSPARSE.CuSparseVector(x)
    gpu!(x::S) where S<:𝕊² = return x = CUDA.CUSPARSE.CuSparseMatrixCSR(x)

    cpu!(x::S) where S<:𝕊¹ = return x = SparseArrays.SparseVector(x)
    cpu!(x::S) where S<:𝕊² = return x = SparseArrays.SparseMatrixCSR(x)
    # undo sparsity while preserving memory location
    dense!(x::S) where S<:CUDA.CUSPARSE.AbstractCuSparseArray = return x = CUDA.CuArray(x)
    # Handle non-functional gpu!s or no gpus.
elseif !has_cuda_gpu() & !has_cuda()
    gpu! = identity        
    cpu! = identity
    @warn "CUDA gpu! found but no has_cuda()! Printing CUDA info. CPU mode."
    CUDA.versioninfo()
elseif !has_cuda_gpu() & has_cuda()
    gpu! = identity        
    cpu! = identity
    @warn "CUDA found but has_cuda_gpu!() failed! Printing CUDA info. CPU mode."
    CUDA.versioninfo()
else
    gpu! = identity        
    cpu! = identity
    @info "No CUDA found. CPU only mode."
end

# gpu -> cpu interface
"""cpu!(s::AbstractSpectrum; recurse = false) -> s
Transfer the top level of a spectrum struct to the CPU & RAM. If recurse = true, also transfer any 
contained spectra such as raw or processed Bruker files.
"""
cpu!(w::WrappedSpectrum; recurse = false) = begin
    ft = cpu!(w.ft)
    fs = cpu!(w.fs)
    recurse && cpu!(w.src; recurse = true)
    w = WrappedSpectrum(w.name, ft, fs, w.src, w.expno)
    return w
end
cpu!(b::BrukerSpectrum; recurse = false) = begin
    ft = cpu!(b.fid)
    b = BrukerSpectrum(ft, b.acqu, b.procs, b.default_proc, b.name, b.expno)
    if recurse 
        [cpu!(b.procs[p]) for p in keys(b.procs)]
    end
    return b
end
cpu!(p::ProcessedSpectrum; recurse) = begin
    im_ft = cpu!(p.im_ft)
    re_ft = cpu!(p.re_ft)
    p = ProcessedSpectrum(re_ft, im_ft, p.params, p.intrng, p.procno, p.title)
    return p
end

# cpu -> gpu interface
"""gpu!(s::AbstractSpectrum; recurse = false) -> s
Transfer the top level of a spectrum struct to the GPU & VRAM. If recurse = true, also transfer any 
contained spectra.
"""
gpu!(w::WrappedSpectrum; recurse = false) = begin
    ft = gpu!(w.ft)
    fs = gpu!(w.fs)
    recurse && gpu!!(w.src; recurse = true)
    w = WrappedSpectrum(w.name, ft, fs, w.src, w.expno)
    return w
end
gpu!(b::BrukerSpectrum; recurse = false) = begin
    ft = gpu!(b.fid)
    b = BrukerSpectrum(ft, b.acqu, b.procs, b.default_proc, b.name, b.expno)
    if recurse 
        [gpu!!(b.procs[p]) for p in keys(b.procs)]
    end
    return b
end
gpu!(p::ProcessedSpectrum; recurse) = begin
    im_ft = gpu!(p.im_ft)
    re_ft = gpu!(p.re_ft)
    p = ProcessedSpectrum(re_ft, im_ft, p.params, p.intrng, p.procno, p.title)
    return p
end
export gpu!, cpu!, dense!
