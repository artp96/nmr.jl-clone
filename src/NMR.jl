module NMR

using DataFrames 
using GLMakie
using LaTeXStrings
using Interpolations
using LinearAlgebra: norm, normalize, dot
using Statistics
using RecipesBase
using FFTW, AbstractFFTs


import Base: show, dump, /, +, *, -, copy!, in
# include structs first
include("spectrum_structs.jl")

include("utils.jl")
include("baseline.jl")
include("bruker/read.jl")
include("bruker/write.jl")
include("composition.jl")
include("decomposition.jl")
include("integration.jl")
include("interpolation.jl")
# old Plots.jl functions
#include("plotting.jl")
include("show.jl")

# AP's things
include("makie_funcs.jl")
include("slice_size.jl")
include("bruker/popt.jl")
include("eachfibre.jl")

try
    if debug
        return debug
    end
catch 
    global debug = false
end

export Spectrum, ProcessedSpectrum, dump, # constructors, I/O
       plot, plot!, # plotting
       lsq_analyze, candidates, decompose, # decomposition
       ppmtoindex, hztoindex, ppmtohz, # unit conversion
       baseline_correct!, # processing
       limits, chemical_shifts, extract, copy!, # utility functions
       interpolate, resample, # interpolation
       integrate, # integration
       intrng, intrng_data, intrng_indices, intrng_shifts, # integration regions
       debug # debug variable
end
