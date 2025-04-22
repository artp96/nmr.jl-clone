# Bruker I/O functions

"""     read_bruker_binary(fname)
Bruker FID, and processed frequency-domain data files (fid, 1r, 1i) contain a flat
binary representation of their data points, each of which is an Int32."""
function read_bruker_binary(fname)
    reinterpret(Int32, open(read, fname))
end

function parse_float_list(m)
    lines = split(m)
    # The first line of a match is "(0..n)"; we don't need this.
    [parse(Float64, s) for s in split(m)[2:end]]
end

function parse_or(::Type{T}, s, default::T) where T
    v = tryparse(T, s)
    v === nothing ? default : v
end

const filters = [ ( Set(["SW", "SW_h", 
    "O1", "O2", "O3", 
    "SFO1", "SFO2", "SFO3", "SF", 
    "BF1", "BF2", "BF3"]),
              s -> parse(Float64, s) ),
            ( Set(["TD", "NS", "DS", "SI"]),
              s -> parse(Int, s) ),
            ( Set(["D", "P", "GPX", "GPY", "GPZ"]),
              parse_float_list),
            ( Set(["PULPROG"]),
              x -> strip(x)[2:end-1] ),
            ( Set(["AUTOPOS"]),
              x -> parse_or(Int, strip(x)[2:end-1], 0) )
]

function read_params(file)
    contents = read(file, String)
    matches = eachmatch(r"##\$?(.*?)=\s?([^#]+)"s, contents)
    res = Dict(string(m.captures[1]) => parse_param(m.captures[1], m.captures[2]) for m in matches)

    # a few little tweaks
    if "O1" in keys(res) && !("O1P" in keys(res))
        res["O1P"] = res["O1"] / res["SFO1"]
    end
    res
end

function read_intrng(file)
    lines = try
        readlines(file)
    catch
        return missing
    end
    if length(lines) > 2 && strip(lines[1])[1] == 'A'
        [tuple(map(s -> parse(Float64, s), split(line)[1:2])...) for line in lines[3:end]]
    elseif length(lines) > 1 && strip(lines[1])[1] == 'P'
        [tuple(map(s -> parse(Float64, s), split(line))...) for line in lines[2:end]]
    else
        return missing
    end
end

function parse_param(param, val)
    for (names, fun) in filters
        if param in names
            return fun(val)
        end
    end
    return strip(string(val))
end

function ProcessedSpectrum(path :: AbstractString, procno :: Int)
    re_ft = float(read_bruker_binary(joinpath(path, "1r")))
    # mc spectra don't have an imag part
    im_ft = zeros(length(re_ft))
    try
        im_ft = float(read_bruker_binary(joinpath(path, "1i")))
    catch e
        @info "$path appears to be an MC or PS spectrum."
    end
    params = read_params(joinpath(path, "proc"))
    title = read(joinpath(path, "title"), String)
    intrng = read_intrng(joinpath(path, "intrng"))
    @debug println(intrng)
    return NMR.ProcessedSpectrum(re_ft, im_ft, params, intrng, procno, title)
end

ProcessedSpectrum(path::AbstractString) = ProcessedSpectrum(path, parse(Int, basename(path)))

"""
    Spectrum("/data/path") -> s :: {Spectrum <: AbstractSpectrum}
----------------------------------------------------------------------------------
Outermost constructor to import a spectrum. Prompts for a procno to set as the 
default procno. Returns a concrete "Spectrum" container type, where aqpars and 
procpars can be accessed using dict indexing, e.g.
    - s["TD"] -> Number of points in the FID, an acqupar.
    - s["LB"] -> Line broadening applied, Hz, a procpar.
"""
Spectrum(path :: AbstractString, procnos :: AbstractArray{Int}, default_proc :: Int) = begin
    # below, changed joinpath to omit "fid", can't find this anywhere ?
    fid = float(read_bruker_binary(joinpath(path, "fid")))
    # fid = float(read_bruker_binary(path))
    acqu = read_params(joinpath(path, "acqu"))
    name = basename(dirname(path))
    expno = basename(path) 

    # Julia's basename function will not work if the path ends in a trailing slash
    if isempty(expno)
        @warn "Julia's main.basename() returns empty if path ends in \"/\". 
              \n Check the path variable ends in the expno, an integer."
        expno = basename(path[1:end-1])
        expno = parse(Int,expno) 
    else        
        expno = parse(Int,expno) 
    end

    # Needs to be instantiated correctly
    procs = Dict{Int, ProcessedSpectrum{Float64, Vector{Float64}}}()

    for procno in procnos
        proc_path = joinpath(path, "pdata", string(procno))
        procs[procno] = ProcessedSpectrum(proc_path, procno)
    end
    Spectrum(fid, acqu, procs, default_proc, name, expno)
end

Spectrum(path :: AbstractString, procno :: Int) = Spectrum(path, [procno], procno)


"""
    Spectrum("/data/path") -> s :: {Spectrum <: AbstractSpectrum}
----------------------------------------------------------------------------------
Outermost constructor to import a spectrum. Prompts for a procno to set as the 
default procno. Returns a concrete "Spectrum" container type, where aqpars and 
procpars can be accessed using dict indexing, e.g.
    - s["TD"] -> Number of points in the FID, an acqupar.
    - s["LB"] -> Line broadening applied, Hz, a procpar.
"""
function Spectrum(path :: AbstractString; interactive = true)
    procpath = joinpath(path, "pdata")
    if interactive
        println.(parse.(Int, readdir(procpath)))
        println("Choose a procno")
        procno = parse.(Int, readline())
    else
        procno = parse.(Int, readdir(procpath))
    end
    Spectrum(path,  procno)
end
Spectrum(path :: AbstractString, procnos :: AbstractArray{Int}) = Spectrum(path, procnos, minimum(procnos))

function multiimport(fpath) 
    N = readdir(fpath)
    data = Vector{Union{Missing,Spectrum}}(undef, length(N))
    fill!(data,missing)
    for n in N
        try       
            data[parse(Int,n)] = Spectrum(joinpath(fpath, n); interactive = false)
        catch e
            println("Could not parse expno $n:\n$e.")
        end
    end
    return data = filter(!ismissing, data)
end
export read_bruker_binary, multiimport
