
for op in (:*, :/)
    @eval begin
        function $op(a::T, s::BrukerSpectrum) where T <: Number
            t = deepcopy(s)
            broadcast!($op, t.fid, t.fid, a)
            for (_, p) in t.procs
                broadcast!($op, p.re_ft, p.re_ft, a)
                broadcast!($op, p.im_ft, p.im_ft, a)
            end
            t
        end
    end
end

for op in (:+, :-, :*)
    @eval begin
        function $op(args::Vararg{BrukerSpectrum, N}) where N
            shifts = union_shifts(collect(args))
            vals = broadcast($op, (resample(s, shifts) for s in args)...)
            # TODO: make a proper processed BrukerSpectrum
            (shifts, vals)
        end
    end
end

function copy!(src::BrukerSpectrum, dst::BrukerSpectrum, shift_rng::Tuple{Float64,Float64})
    d = resample(src, chemical_shifts(dst))
    rng = ppmtoindex(dst, shift_rng)
    dst[rng] = d[rng]
    dst
end
