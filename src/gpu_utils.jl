
# cpu -> gpu interface
"""cpu!(s::AbstractSpectrum; recurse = false) -> s
Transfer a spectrum struct to the cpu & RAM. If recurse = true, also transfer any 
contained spectra.
"""
cpu!(w::WrappedSpectrum; recurse = false) = begin
    ft = cpu(w.ft)
    fs = cpu(w.fs)
    recurse && cpu!(w.src; recurse = true)
    w = WrappedSpectrum(w.name, ft, fs, w.src, w.expno)
    return w
end
cpu!(b::BrukerSpectrum; recurse = false) = begin
    ft = cpu(b.fid)
    b = BrukerSpectrum(ft, b.acqu, b.procs, b.default_proc, b.name, b.expno)
    if recurse 
        [cpu!(b.procs[p]) for p in keys(b.procs)]
    end
    return b
end
cpu!(p::ProcessedSpectrum; recurse) = begin
    im_ft = cpu(p.im_ft)
    re_ft = cpu(p.re_ft)
    p = ProcessedSpectrum(re_ft, im_ft, p.params, p.intrng, p.procno, p.title)
    return p
end

# gpu -> cpu interface
"""gpu!(s::AbstractSpectrum; recurse = false) -> s
Transfer a spectrum struct to the gpu & VRAM. If recurse = true, also transfer any 
contained spectra.
"""
gpu!(w::WrappedSpectrum; recurse = false) = begin
    ft = gpu(w.ft)
    fs = gpu(w.fs)
    recurse && gpu!(w.src; recurse = true)
    w = WrappedSpectrum(w.name, ft, fs, w.src, w.expno)
    return w
end
gpu!(b::BrukerSpectrum; recurse = false) = begin
    ft = gpu(b.fid)
    b = BrukerSpectrum(ft, b.acqu, b.procs, b.default_proc, b.name, b.expno)
    if recurse 
        [gpu!(b.procs[p]) for p in keys(b.procs)]
    end
    return b
end
gpu!(p::ProcessedSpectrum; recurse) = begin
    im_ft = gpu(p.im_ft)
    re_ft = gpu(p.re_ft)
    p = ProcessedSpectrum(re_ft, im_ft, p.params, p.intrng, p.procno, p.title)
    return p
end
export gpu!, cpu!
