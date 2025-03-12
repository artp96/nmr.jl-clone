# A function to convert a popt array to an n-dimensional spectrum.
# Spectrum(path :: AbstractString) = begin
#     procnos = [parse(Int,n) for n in readdir(joinpath(path, "pdata"))]
#     popt_detect = any(diff(procnos) > 1)
#     popt_detect ? return Popt_Spectrum(path, procnos) : return Spectrum(path, procnos)
# end

Popt_Spectrum(path :: AbstractString, procnos :: AbstractArray{Int}, default_proc :: Int) = begin
    # below, changed joinpath to omit "fid", can't find this anywhere ?
    #fid = float(read_bruker_binary(joinpath(path, "fid")))
    fid = float(read_bruker_binary(path))
    acqu = read_params(joinpath(path, "acqu"))
    name = basename(dirname(path))
    expno = parse(Int, basename(path))
    procs = Dict()
    for procno in procnos
        proc_path = joinpath(path, "pdata", string(procno))
        procs[procno] = ProcessedSpectrum(proc_path, procno)
    end
    Spectrum(fid, acqu, procs, default_proc, name, expno)
end

