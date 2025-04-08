# script to test heatmap script
using NMR
#    pkg needed to resolve dependencies while developing


using BenchmarkTools
"""
    Some useful globals
"""
global debug = false
spath  = joinpath(homedir(), "nmrdata/data/2025/b300b10/250318")
expno = "7"
#spath  = joinpath(homedir(), "nmrdata/data/2025/b500b07/250326")
#expno = "9"

#!TODO:  NMR.Spectrum does not work without 1i/1r for now 130325. Need to test on a zg spectrum. 

#!TODO:  NMR.Spectrum does not work without 1i/1r for now 130325. Need to test on a zg spectrum. 
# Need to test on a zg spectrum

s = Spectrum(joinpath(spath, expno), [1], 1)

poptpath = joinpath(spath, expno)
protpath = joinpath(poptpath, "popt.protocol.999")

# testing inside of parse_popt()
file = protpath; contents = read(file, String)

df = read_PoptProtocol(file);
dims, vars = get_ArrayPoptDims(df)
mat = restructure_array(df, df.var"Maximum point");

using GLMakie
GLMakie.activate!()
#heatmap(mat)
#! TODO need to get the dims and variables from this function instead, so restructure_array() wraps an 
# get_array_structure(df) -> (dims, axis_variables) type setup, then
fig = splatted_heatmaps(df);
#

S = PoptSpectrum(poptpath; UI_enable = true);
using BenchmarkTools
import NMR: trim_spectrum
fibs = eachfibre(S.ft * S.ser)
idxs = trim_spectrum(S, 2.104)
b1 = @benchmark p₁ = auto_ϕ_correct(fibs[1], idxs)
b2 = @benchmark p₂ = auto_ϕ_correct2(fibs[1], idxs)

