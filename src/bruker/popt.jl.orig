

PoptSpectrum(path :: AbstractString, procno :: Int; kwargs...) = Popt_Spectrum(path, [procno], procno; kwargs...)

"""
    PoptSpectrum(path :: AbstractString; UI_enable = true)
First outer constructor for simple user interface mode."""
function PoptSpectrum(path :: AbstractString; UI_enable = true)
    procnos = parse.(Int, readdir(path * "/pdata"))
    if UI_enable
        println("Choose a procno as the target popt array, return for auto.")
        println.(procnos)
        println("---")
        procno = readline()
    end
    if !UI_enable || isempty(procno) 
        procno = maximum(procnos)
    end
    procno isa String ? procno = parse(Int, procno) : procno = procno
    PoptSpectrum(path, procnos, procno; UI_enable = UI_enable)
end

"""
    PoptSpectrum(path, procnos; kwargs...) -> PoptSpectrum
    path :: AbstractString,
    procnos :: AbstractArray{Int}
____________________________________________________________________________
Second outer constructor wrapper.
"""
PoptSpectrum(path :: AbstractString, procnos :: AbstractArray{Int}; kwargs...) = Popt_Spectrum(
    path, procnos, maximum(procnos); kwargs... )

""" 
    PoptSpectrum(path, procnos, procno; poptno, UI_enable) -> PoptSpectrum
    path :: AbstractString, 
    procnos :: AbstractArray{Int},
    procno :: Int;
    poptno = "", UI_enable = true)
____________________________________________________________________________
Final outer constructor to Struct, containing a multidimensional spectrum 
derived from a POPT array.
"""
function PoptSpectrum(path :: AbstractString, procnos :: AbstractArray{Int}, procno :: Int; poptno = "", UI_enable = true) 
    !isempty(poptno) && prepend!(poptno, ".")
    # below, changed joinpath to omit "fid", can't find this anywhere ?
    #fid = float(read_bruker_binary(path))
    acqu = read_params(joinpath(path, "acqu"))
    protocol = read_PoptProtocol(joinpath(path, "popt.protocol" * poptno))
    ser = fetch_serfile(path, poptno, contents; UI_enable = UI_enable)
    name = basename(dirname(path))
    expno = parse(Int, basename(path))  
    # Popt spectra arrays need only 1 procno and 1 serfile
    proc_path = joinpath(path, "pdata", string(procno))
    proc = ProcessedSpectrum(proc_path, procno)
    return PoptSpectrum(acqu, procno, expno, name, path, protocol, ser, proc)
end

"""
    fetch_serfile(path, poptno, contents; UI_enable = false) -> Vector{f64}
____________________________________________________________________________
Function to get the ser matched to the given popt.protocol.n
# TODO Not working properly with UI_enable off, need to pass right match target.
"""
function fetch_serfile(path, poptno, contents; UI_enable = false)
    path = splitpath(path)
    path = joinpath(path[1:end-1])    
    #prompt for the correct expno
    if isempty(poptno) && UI_enable
        println("Choose an expno as the popt serfile container.")
        println.(parse.(Int, readdir(path)))
        println("---")
        expno = readline()
        ser_addr = joinpath(path, expno, "ser")
        ser = read_bruker_binary(ser_addr) |> float
    elseif !isempty(poptno)
        expno = poptno
        ser_addr = joinpath(path, expno, "ser")
        ser = read_bruker_binary(ser_addr)
    elseif isempty(poptno) & !UI_enable
        m = match(r"serfile: (.*)\s", contents)
        try
            expno = splitpath(m[1])[end]
            ser_addr = joinpath(path, expno, "ser")
            ser = read_bruker_binary(ser_addr) |> float
        catch err
            @warn "No expno found for serfile."
            ser = Vector{typeof(1.)}([]) # not an unusual case, handle gracefully
        end
    end
    return ser
end

"""
    parse_PoptProtocol(file; run = 1) -> headers, cols
____________________________________________________________________________
Function to parse the popt protocol into a set of headers and only 
interesting columns.
Wrapped by read_PoptProtocol to give a DataFrame
"""
function parse_PoptProtocol(file; run = 1)
    contents = read(file, String)
    debug && println(contents)
    # Extract headers 
    headerline = eachmatch(r"MOD=\s*\S+\n\n(Experiment[A-Za-z\d\s]+ Integral)\n\s*\d"s,
                           contents)
    length(collect(headerline)) < 1 && @error "No matches found."
    headerline = collect(headerline)[run]

    headers = filter!(!isempty, split(headerline.captures[1], r"\s\s"))
    # remove leading/trailing whitespace
    headers = strip.(headers)

    #headers = hcat([h.captures for h in headers]...)
    # col_regex = build_data_regex(headers)

    # get all the data between the headers and the footer 
    columns = eachmatch(r"Integral\n([\s\d\-\.\n]+)\npoptau"s, contents)
    columns = collect(columns)[run]
    # separate by rowin
    columns = split(columns.captures[1], "\n")
    # split between spaces, drop empty strings:
    columns = split.(columns, r"\s+")
    for c in columns 
        filter!(!isempty, c)
    end
    # Splat as a matrix of substrings into numbers
    columns = parse.(Float64, hcat(columns...))
    return (headers, columns)
end

"""
    read_PoptProtocol(file) -> DataFrame
____________________________________________________________________________
Function to construct a table of popt.protocol data 
"""
function read_PoptProtocol(file)
    h, c = parse_PoptProtocol(file)
    df = DataFrame(transpose(c), vec(h))
    # for Type safety
    df.Experiment = convert.(Int64, df.Experiment)
    return df
end

"""
    build_data_regex(headers <: Base.RegexMatchIterator{String}) -> Regex
____________________________________________________________________________
Form a dynamic regex to split the columns of data - redundant?
"""
function build_data_regex(headers :: T) where T <: Base.RegexMatchIterator{String}
    #Can't natively get the length of the Iterator type
    count = 0
    for _ in headers
        count += 1
    end
    regex_str = "^" * join(repeat(["(\\S+)"], count), "\\s+") * "\$"
    return Regex(regex_str)
end

"""
    get_ArrayPoptDims(df <: AbstractDataFrame) -> dims, axes

____________________________________________________________________________
A function to scan the dependent variable columns and deduce the array 
structure, yielding the dimensions of the Popt experiment and the x,y,z 
axes a NamedTuple of (:Variable, [Values]) pairs.
"""
function get_ArrayPoptDims(df :: T) where T <: AbstractDataFrame
    boring_columns = ("Experiment", "Maximum point", "Minimum point", "Integral")
    # find the columns with dependent variables, and filter the others out
    h = names(df)
    filter!(x -> !(x ∈ boring_columns), h)

    # Want the number of independent variables, excluding N, min, max, ∫dx cols.
    n = length(h)
    axes = df[:, 2:n+1]
    # Find the unique combinations of variables
    vals = unique.(eachcol(axes))
    # must be a Tuple
    dims = length.(vals)
    debug && println(dims)
    indices = NamedTuple(zip(Symbol.(h), vals))
    
    # popt always does the last variable first, which julia sees as the reverse of column-
    # major indexing, so the returned variables should be reversed
    reverse!(dims)
    dims = tuple.(dims...)
<<<<<<< HEAD

    return dims, axes
=======
    # reversing a named tuple only reverses the key order, which is appropriate 
    idxs = reverse(indices)
    return dims, idxs
>>>>>>> 31b7c2e9a81f2870f2854a711fd8d9f60d60ec55
end

"""
    restructure_array(df::DataFrame, vec) -> reshape(vec)
____________________________________________________________________________
A function to take an array of integrals or intensities from experiments
and arrange them according to the popt protocol used.
Should work to take f(df, df.int), or f(df, vec) for some custom vector -
the purpose of this function is to order indices.
"""
function restructure_array(df::DataFrame, vec)
    #get the new array structure
    dims, vars = get_ArrayPoptDims(df)
    # return the reshaped array
    return reshape(vec, dims)
end

export read_PoptProtocol, parse_PoptProtocol, restructure_array, PoptSpectrum, get_ArrayPoptDims, fetch_serfile
