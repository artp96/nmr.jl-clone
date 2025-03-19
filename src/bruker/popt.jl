# A function to convert a popt array to an n-dimensional spectrum.
# Spectrum(path :: AbstractString) = begin
#     procnos = [parse(Int,n) for n in readdir(joinpath(path, "pdata"))]
#     popt_detect = any(diff(procnos) > 1)
#     popt_detect ? return Popt_Spectrum(path, procnos) : return Spectrum(path, procnos)
# end
using DataFrames 

""" Popt_Spectrum()
    Struct to containing a multidimensional spectrum derived from a POPT 
    array.
"""
Popt_Spectrum(path :: AbstractString, procnos :: AbstractArray{Int}, default_proc :: Int) = begin
    # below, changed joinpath to omit "fid", can't find this anywhere ?
    #fid = float(read_bruker_binary(joinpath(path, "fid")))
    fid = float(read_bruker_binary(path))
    acqu = read_params(joinpath(path, "acqu"))
    serfile = fetch_serfile(path)
    protocol = read_PoptProtocol(joinpath(path, "popt.protocol"))
    name = basename(dirname(path))
    expno = parse(Int, basename(path))  
    procs = Dict()
    popt_protocol = read_PoptProtocol(joinpath(path, "popt.protocol"))
    for procno in procnos
        proc_path = joinpath(path, "pdata", string(procno))
        procs[procno] = ProcessedSpectrum(proc_path, procno)
    end
    Spectrum(fid, acqu, procs, default_proc, name, expno)
end

function fetch_serfile(file)
    contents = read(file, String)
    try  
    # This regex captures the end of this line to the new line.
      serfile_addr = match(r"Target directory for serfile: (.*?)\n", contents).captures[1]
    catch
      # can handle with !ismissing later
      serfile_addr = missing
    end     
    return serfile_addr
end

"""
    Function to parse the popt protocol into a set of headers and columns
"""
function parse_PoptProtocol(file; run = 1, debug1 = true)
    contents = read(file, String)
    debug1 && println(contents)
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
    Form a dynamic regex to split the columns of data
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
    A function to scan the dependent variable columns and deduce the array structure
"""
function get_ArrayPoptDims(df :: T) where T <: AbstractDataFrame
    boring_columns = ("Experiment", "Maximum point", "Minimum point", "Integral")
    # find the columns with dependent variables
    h = names(df)
    filter!(x -> !(x ∈ boring_columns), h)

end

"""
    A function ake an array of integrals or intensities from experiments
    and arrange them according to the popt protocol used.
    Should work to take f(df, df.int), or f(df, vec) for some custom vector -
    the purpose of this function is to order indices.
"""
function restructure_array(popt_Table::DataFrame, vec)
    h = names(popt_Table)
    # Want the number of independent variables, excluding N, min, max, ∫dx cols.
    n = length(h) - 4
    vars = df[:, 2:n+1]
    # Find the unique combinations of variables
    xyz = unique.(eachcol(vars))
    # must be a Tuple
    dims = length.(xyz)
    reverse!(dims)
    println(dims)
    dims = tuple.(dims...)
    # return the reshaped array
    return reshape(vec, dims)
end



export read_PoptProtocol, parse_PoptProtocol, restructure_array
