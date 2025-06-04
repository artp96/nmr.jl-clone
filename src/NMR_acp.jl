module NMR_acp

"""γₚ, ¹H gyromagnetic ratio in MHz T⁻¹. """
const γₚ = 42.577478461 
const dpath = joinpath(homedir(), "nmrdata/data/2025")

using DataFrames 
using DSP # Digital Signal processing
using ProgressBars
using DataFrames 
using GLMakie; GLMakie.activate!()
using LaTeXStrings
using Interpolations
using LinearAlgebra: norm, normalize, dot
using Statistics
using RecipesBase
using FFTW
using FastTransforms
using FFTW
using ForwardDiff
using GLMakie; GLMakie.activate!()
using Interpolations
using LaTeXStrings
using LinearAlgebra: norm, normalize, dot
using Optimization
using Optim
using ProgressBars
using Statistics
using StatsBase
using RecipesBase

using CUDA
if has_cuda() & has_cuda_gpu()
    @info "CUDA found. GPU mode."
    gpu = CuArray
    cpu = Array
elseif has_cuda_gpu() & !has_cuda()
    gpu = identity        
    cpu = identity
    @warn "CUDA gpu found but no has_cuda()! Printing CUDA info."
    CUDA.versioninfo()
elseif !has_cuda_gpu() & has_cuda()
    gpu = identity        
    cpu = identity
    @warn "CUDA found but has_cuda_gpu() failed! Printing CUDA info."
    CUDA.versioninfo()
else
    @info "No CUDA found. CPU mode."
end

export gpu, cpu

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
# using Plots
#include("plotting.jl")
include("show.jl")
include("gpu_utils.jl")

# AP's things
include("makie_funcs.jl")
include("slice_size.jl")
include("bruker/popt.jl")
include("bruker/fid_import.jl")
include("eachfibre.jl")
include("processing/phase_correction.jl")
include("processing/window_functions.jl")
include("processing/transforms.jl")
include("processing/zerofill.jl")

try
    if debug
        return debug
    end
catch 
    global debug = false
end

export BrukerSpectrum, ProcessedSpectrum, dump, # constructors, I/O
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
