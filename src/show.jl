
function show(io::IO, s::BrukerSpectrum)
    println(io, "==== $(s.name) ($(s.expno)) ====")
    println(io, "FID size: $(length(s.fid))")
    println(io, "# processings: $(length(s.procs))") 
    println(io, "Default proc. no.: $(s.default_proc)")
    out = PipeBuffer()
    for k in sort(collect(keys(s.procs)))
        show(out, s.procs[k])
        println(out, "")
    end
    for l in eachline(out)
        println(io, "    $l")
    end
end

function show(io::IO, p::ProcessedSpectrum)
    println(io, "---- Proc. no. $(p.procno) ----")
    println(io, strip(p.title))
    println(io, "Processed size: $(p["SI"]) points")
end

function show(io::IO, w::FracSpectrum)
    println(io, "f(s): $(typeof(w.fs)), size: $(size(w.fs)...)")
    println(io, "f(t): $(typeof(w.ft)), size: $(size(w.fs)...)")
    println(io, "\n vvv  W.src  vvv")
    println(io, "==== $(w.src.name) ($(w.src.expno)) ====")
    out = PipeBuffer()
    for k in sort(collect(keys(w.src.procs)))
        show(out, w.src.procs[k])
        println(out, "")
    end
    for l in eachline(out)
        println(io, "    $l")
    end
end
