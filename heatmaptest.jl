# script to test heatmap script
#
include("src/NMR.jl")
include("src/bruker/popt.jl")

nmrdata = "~/nmrdata/data" 
spath  = "/home/w25612ap/nmrdata/data/AP/2025/b300b10/250318"

#!TODO:  NMR.Spectrum does not work without 1i/1r for now 130325. Need to test on a zg spectrum. 
# Need to test on a zg spectrum


S = NMR.Spectrum(spath * "/7", [1], 1)

poptpath = spath * "/7";
protpath = poptpath * "/popt.protocol.999";

# testing inside of parse_popt()
file = protpath; contents = read(file, String)

df = read_PoptProtocol(file)

mat = restructure_array(df, df.var"Maximum point")

using GLMakie
GLMakie.activate!()
#heatmap(mat)

fig = splatted_heatmaps(mat)
