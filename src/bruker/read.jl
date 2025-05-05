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

#TODO: Make const
#const
filters = [ ( Set(["SW", "SW_h", 
    "O1", "O2", "O3", 
    "SFO1", "SFO2", "SFO3", "SF", 
    "BF1", "BF2", "BF3", 
    "PHC0", "PHC1", "LB", 
    # GRPDLY / DSP Variables.
    "GRPDLY", "DSPFVS", "DECIM"]),
              s -> parse(Float64, s) ),
            ( Set(["TD","TD0", "NS", "DS", "SI", "NC", "NC_proc",
                "FnMODE",]),
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
    
    params = read_params(joinpath(path, "proc"))
    
    if "1r" in readdir(path)
        (r,i) = ("1r", "1i")
    elseif "2rr" in readdir(path)
        (r,i) = ("2rr", "2ii")
    else
        @error "No 1r or 2rr file found for $procno \n (at $path)."
    end
        re_ft = read_bruker_binary(joinpath(path, r)) .* 2^params["NC_proc"]

    im_ft = zeros(size(re_ft))
    try 
        im_ft = read_bruker_binary(joinpath(path, i)) .* 2^params["NC_proc"]
    catch e
        @info "$path has no imaginary part."
    end
    scale_correction = params
    title = read(joinpath(path, "title"), String)
    intrng = read_intrng(joinpath(path, "intrng"))
    @debug println(intrng)
    return NMR.ProcessedSpectrum(float(re_ft), float(im_ft), params, intrng, procno, title)
end

ProcessedSpectrum(path::AbstractString) = ProcessedSpectrum(path, parse(Int, basename(path)))

"""
   BrukerSpectrum("/data/path") -> s :: {BrukerSpectrum <: AbstractSpectrum}
----------------------------------------------------------------------------------
Outermost constructor to import a spectrum. Prompts for a procno to set as the 
default procno. Returns a concrete "Spectrum" container type, where aqpars and 
procpars can be accessed using dict indexing, e.g.
    - s["TD"] -> Number of points in the FID, an acqupar.
    - s["LB"] -> Line broadening applied, Hz, a procpar.
"""
BrukerSpectrum(path :: AbstractString, procnos :: AbstractArray{Int}, default_proc :: Int) = begin
    # below, changed joinpath to omit "fid", can't find this anywhere ?
    # fid = float(read_bruker_binary(path))
    acqu = read_params(joinpath(path, "acqus"))
    
    fid = zeros(acqu["TD"])
    if "fid" in readdir(path)
        fid = float(read_bruker_binary(joinpath(path, "fid"))) 
    elseif "ser" in readdir(path)
        fid = float(read_bruker_binary(joinpath(path, "ser"))) 
    else
        @error "No 'fid' or 'ser' found in $path."
    end
    fid .*= 2.0^acqu["NC"]
     
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
   BrukerSpectrum(fid, acqu, procs, default_proc, name, expno)
end

BrukerSpectrum(path :: AbstractString, procno :: Int) = BrukerSpectrum(path, [procno], procno)


"""
   BrukerSpectrum("/data/path") -> s :: {BrukerSpectrum <: AbstractSpectrum}
----------------------------------------------------------------------------------
Outermost constructor to import a spectrum. Prompts for a procno to set as the 
default procno. Returns a concrete "Spectrum" container type, where aqpars and 
procpars can be accessed using dict indexing, e.g.
    - s["TD"] -> Number of points in the FID, an acqupar.
    - s["LB"] -> Line broadening applied, Hz, a procpar.
"""
function BrukerSpectrum(path :: AbstractString; interactive = true)
    procpath = joinpath(path, "pdata")
    procnos = parse.(Int, readdir(procpath))
    if !interactive || length(procnos) == 1
        default_proc = minimum(procnos)
    else
        println.(procnos)
        println("Choose a procno")
        default_proc = parse.(Int, readline())
    end
    BrukerSpectrum(path, procnos, default_proc)
end

BrukerSpectrum(path :: AbstractString, procnos :: AbstractArray{Int}) = BrukerSpectrum(path, procnos, minimum(procnos))

function multiimport(fpath) 
    N = readdir(fpath)
    data = Vector{Union{Missing, BrukerSpectrum}}(undef, length(N))
    fill!(data,missing)
    for n in N
        try       
            data[parse(Int,n)] = BrukerSpectrum(joinpath(fpath, n); interactive = false)
        catch e
            println("Could not parse expno $n:\n$e.")
        end
    end
    # ensures element type stability
    return skipmissing(data) |> collect
end
export read_bruker_binary, multiimport
