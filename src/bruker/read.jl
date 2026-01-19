# Bruker I/O functions

"""     read_bruker_binary(fname)
Bruker FID, and processed frequency-domain data files (fid, 1r, 1i) contain a flat
binary representation of their data points, each of which is an Int32."""
function read_bruker_binary(fname)
    reinterpret(Int32, open(read, fname))
end

function clean_match(s)
    strip(replace(s, r"[()\"]" => ""))
end
function parse_list(T::Type, m)
    lines = filter(!isempty,clean_match.(split(m)))
    # The first line of a match is "(0..n)"; we don't need this.
    [parse(T, s) for s in lines[2:end]]
end

function parse_or(::Type{T}, s, default::T) where T
    v = tryparse(T, s)
    v === nothing ? default : v
end

const filters = [ 
    (
        Set(["SW", "SW_h", 
        "O1", "O2", "O3", "O4",
        "SFO1", "SFO2", "SFO3", "SFO4", 
        "SF",
        "BF1", "BF2", "BF3", 
        "PHC0", "PHC1", "LB", 
        # GRPDLY / DSP Variables.
        "GRPDLY", "DSPFVS", "DECIM",
        "DWELL", "DE", "TE", "NusAMOUNT"]),
        s -> parse(Float64, clean_match(s)) 
    ),

    # Int64eger parameters, such as number of scans, and some enum keys
    (
        Set(["TD","TD0", "NS", "DS", "SI", "NC", "NC_proc",
        "FnMODE", "SEOUT"]),
        s -> parse(Int, clean_match(s) ) 
    ),
    # Lists of floating point numbers, such as pulse-specific parameters
    (
        Set(["D", "P", "GPX", "GPY", "GPZ", "CNST", "SPW", "SPDB", "PLW",
    "PLDB", "SPOFFS", "SPOAL"]),
        s -> parse_list(Float64, s) 
    ),
    ( 
        Set(["TD_INDIRECT"]),
        s -> parse_list(Int, s) 
    ),
    # Pulse program name
    ( 
        Set(["PULPROG"]),
        x -> strip(x)[2:end-1] 
    ),
    ( 
        Set(["AUTOPOS"]),
        x -> parse_or(Int64, strip(x)[2:end-1], 0) 
    ),
    ( 
        Set(["SPNAM", "GPNAM"]),
        x -> parse_or(Int64, strip(x)[2:end-1], 0) 
    )
]

function read_params(file)
    contents = read(file, String)
    res = Dict{String, Any}()
    matches = eachmatch(r"##\$?(.*?)=\s?([^#]+)"s, contents)
    for m in matches
        parsed = parse_param(m.captures[1], m.captures[2])
        res[string(m.captures[1])] = parsed
    end
    # a few little tweaks
    if "O1" in keys(res) && !("O1P" in keys(res))
        res["O1P"] = res["O1"] / res["SFO1"]
    end
    return res
    # return an empty dict if the file isn't found.
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

function parse_param(param, val, filters = filters)
    for (names, fun) in filters
         param in names && return fun(val)
    end
    return strip(string(val))
end

function ProcessedSpectrum(path :: AbstractString, procno :: Int; no_proc_data = true)
    
    params = read_params(joinpath(path, "procs"))
    
    re_ft = im_ft = Float64[]
    # if skipping import of procnos
    if no_proc_data
        @info "Skipping import of procno spectrum."
    else
        if "1r" in readdir(path)
            (r,i) = ("1r", "1i")
            _err = ""; # empty string
        elseif "2rr" in readdir(path)
            (r,i) = ("2rr", "2ii")
            _err = ""
        else
            _err = "No 1r or 2rr file found for $procno \n (at $path)."
        end
        try 
            re_ft = read_bruker_binary(joinpath(path, r)) .* 2^params["NC_proc"]
        catch e
            @info "Procno $path not readable. \n $_err"
        end
        try 
            im_ft = read_bruker_binary(joinpath(path, i)) .* 2^params["NC_proc"]
        catch e
            @info "$path has no imaginary part."
        end
    end
    scale_correction = params
    title = read(joinpath(path, "title"), String)
    intrng = read_intrng(joinpath(path, "intrng"))
    @debug println(intrng)
    return ProcessedSpectrum(float(re_ft), float(im_ft), params, intrng, procno, title)
end

ProcessedSpectrum(path::AbstractString; kwargs...) = ProcessedSpectrum(path, parse(Int, basename(path)); kwargs...)

"""
merge_acqus([dict1, dict2, ...]) -> dict
----------------------------------------------------------------------------------
Function to merge acquisition parameter dicts. Most of the useful stuff is in 
acqu1, including TD_INDIRECT, but NUS parameters are stored in each acquNs file.
Appends a dimension-specific number to keys in acqu2, acqu3, etc. but not to keys 
in acqu1. 
    acqus["TD"] -> acqus["TD"]
    acqu2s["TD"] -> acqus["TD2"]
"""
function merge_acqus!(acqu, acquN::Vector{D}) where D<:Dict
    for (n, d) in enumerate(acquN)
        for k in keys(d)
            if k ∉ keys(acqu) 
                acqu[k] = d[k]
            else
                newkey = join([k string(n)])
                acqu[newkey] = d[k]
            end
        end
    end
    return acqu
end
merge_acqus(acquN::Vector{D}) where D<:Dict = merge_acqus!(similar(acquN[1]), acquN[1:end])

"""
    Compute the dimension of the data for reshape(). The TD2 dimension is implicitly zero filled 
    once due to serial acquisition of re,im,re data in Bruker format.
"""
get_ser_dims(D::Dict)::Tuple = (2D["TD2"], filter(!iszero, D["TD_INDIRECT"])...)

"""
   BrukerSpectrum("/data/path") -> s :: {BrukerSpectrum <: AbstractSpectrum}
----------------------------------------------------------------------------------
Outermost constructor to import a spectrum. Prompts for a procno to set as the 
default procno. Returns a concrete "Spectrum" container type, where aqpars and 
procpars can be accessed using dict indexing, e.g.
    - s["TD"] -> Number of points in the FID, an acqupar.
    - s["LB"] -> Line broadening applied, Hz, a procpar.
"""
BrukerSpectrum(path :: AbstractString, procnos :: AbstractArray{Int}, default_proc :: Int; no_proc_data = true) = begin
    # below, changed joinpath to omit "fid", can't find this anywhere ?
    # fid = float(read_bruker_binary(path))
    acqu_files = filter(!isnothing, match.(r"acq.+s", readdir(path)))
    # get higher-dimensional acqus and merge them into a single acqus
    acqu = Dict{String, Any}()
    merge_acqus!(acqu, [read_params(joinpath(path, aq.match)) for aq ∈ acqu_files])
    
    if "fid" in readdir(path)
        TD = acqu["TD"]
        fid = zeros(TD) 
        fid = float(read_bruker_binary(joinpath(path, "fid"))) |> gpu!
        fid = split_fid(fid) |> cpu!
        # fid should be zero_filled immediately to get original dimension.
        fid = zero_fill(fid, 2)
    elseif "ser" in readdir(path)
        fid = float(read_bruker_binary(joinpath(path, "ser"))) |> gpu!
        fid = reshape(fid, get_ser_dims(acqu))
        @show typeof(fid)
        fid = split_fid(fid) |> cpu!
        @show typeof(fid)
        fid = zero_fill(fid, 2)
        @show typeof(fid)
    else
        @error "No 'fid' or 'ser' found in $path."
    end
    fid *= 2.0^acqu["NC"]
     
    name = basename(dirname(path))
    expno = basename(path) 

    # Julia's basename function will not work if the path ends in a trailing slash
    if isempty(expno)
        @warn "Julia's main.basename() returns empty if path ends in \"/\". 
              \n Check the path variable ends in the expno, an integer character."
        expno = basename(path[1:end-1])
        expno = parse(Int, expno) 
    else        
        expno = parse(Int, expno) 
    end

    # Needs to be instantiated correctly
    procs = Dict{Int, ProcessedSpectrum{Float64, Vector{Float64}}}()
    
    for procno in procnos
        proc_path = joinpath(path, "pdata", string(procno))
        procs[procno] = ProcessedSpectrum(proc_path, procno; no_proc_data = no_proc_data)
    end
   BrukerSpectrum(fid, acqu, procs, default_proc, name, expno)
end

BrukerSpectrum(path :: AbstractString, procno :: Int; kwargs...) = BrukerSpectrum(path, [procno], procno; kwargs...)


"""
   BrukerSpectrum("/data/path") -> s :: {BrukerSpectrum <: AbstractSpectrum}
-------------------------------------------------------------------------------------------
Outermost constructor to import a spectrum. Prompts for a procno to set as the 
default procno. Returns a concrete "Spectrum" container type, where aqpars and 
procpars can be accessed using dict indexing, e.g.
    - s["TD"] -> Number of points in the FID, an acqupar.
    - s["LB"] -> Line broadening applied, Hz, a procpar.
"""
function BrukerSpectrum(path :: AbstractString; interactive = true, kwargs...)
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

BrukerSpectrum(path :: AbstractString, procnos :: AbstractArray{Int}; kwargs...) = BrukerSpectrum(path, procnos, minimum(procnos); kwargs...)


"""
    multiimport(fpath::AbstractString) 
-------------------------------------------------------------------------------------------
Import an array of Bruker spectra from a folder as FracSpectrum.
- fpath : path/to/folder/of/expnos
- givemissing : Optionally keep the indices of files which 
could not be imported as v[idx] = "missing".
"""
function multiimport(fpath::AbstractString; givemissing = false, kwargs...) 
    N = readdir(fpath)
    data = Vector{Union{Missing, BrukerSpectrum}}(undef, length(N))
    fill!(data,missing)
    Threads.@threads for i in eachindex(data)
        try       
            data[i] = BrukerSpectrum(joinpath(fpath, N[i]; kwargs...); interactive = false)
        catch e
            println("Could not parse expno $(N[i]):\n$e.")
        end
    end
    # ensures element type stability
    if !givemissing 
        filter!(!ismissing, data)
        data = convert(Vector{BrukerSpectrum}, data)
    end
    return data
end 

"""
    multiwrap(fpath::AbstractString) 
-------------------------------------------------------------------------------------------
Import an array of Bruker spectra from a folder as FracSpectrum.
- fpath : path/to/folder/of/expnos
- givemissing : Optionally keep the indices of files which 
could not be imported as v[idx] = "missing".
Calls multiimport under the hood.
"""
multiwrap(V::W; givemissing = false, index = true) where W <: AbstractVector = begin
    wrap = Vector{Union{Missing, FracSpectrum}}(undef, length(V))
    Threads.@threads for i in eachindex(V)
        if !ismissing(V[i])
            wrap[i] = FracSpectrum(V[i]) 
        else
            @info("Expno $i was not imported.")
            wrap[i] = missing
        end
    end
    if !givemissing 
        filter!(!ismissing, wrap) 
        data = convert(Vector{FracSpectrum}, wrap)
    end
    return wrap = index ? Dict([w.expno => w for w in wrap]) : wrap
end
multiwrap(fpath::AbstractString; kwargs...) = multiwrap(multiimport(fpath; kwargs...))

Base.getindex(d::D, r::R) where {D<:Dict{Int,<:AbstractSpectrum}, R<:UnitRange} = [d[i] for i in r]
export read_bruker_binary, multiimport, multiwrap
