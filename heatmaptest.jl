# script to test heatmap script
using NMR
#    pkg needed to resolve dependencies while developing



"""
    Some useful globals
"""
global debug = false
spath  = joinpath(homedir(), "nmrdata/data/2025/b300b10/250318")

#!TODO:  NMR.Spectrum does not work without 1i/1r for now 130325. Need to test on a zg spectrum. 
# Need to test on a zg spectrum


s = Spectrum(spath * "/7", [1], 1)

poptpath = joinpath(spath, "7")
protpath = joinpath(poptpath, "popt.protocol")

# testing inside of parse_popt()
file = protpath; contents = read(file, String)

df = read_PoptProtocol(file);
h = get_ArrayPoptDims(df);
mat = restructure_array(df, df.var"Maximum point");

using GLMakie
GLMakie.activate!()
#heatmap(mat)
#! TODO need to get the dims and variables from this function instead, so restructure_array() wraps an 
# get_array_structure(df) -> (dims, axis_variables) type setup, then
fig = splatted_heatmaps(mat);
#

S = PoptSpectrum(poptpath; UI_enable = true); p.ser

P = auto_ϕ_correct(S)
