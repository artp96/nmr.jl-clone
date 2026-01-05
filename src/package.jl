using BenchmarkTools
using Interpolations
using DataFrames
using SparseArrays
using FastTransforms
using FFTW
using GLMakie; GLMakie.activate!()
using LaTeXStrings
using LinearAlgebra: norm, normalize, dot
using ProgressBars
using RecipesBase
using DSP
using Optimization, ForwardDiff
using Base: debug_color
using OptimizationOptimJL
using CUDA

import Base: show, dump, /, +, *, -, copy!, in

# set some constants at this level

""" 
    Gyromagnetic ratios in MHz T⁻¹. 
"""
const γ = (
    H1 = 42.577478461,# ½ 𝑆, 99.9% abund.
    Pr = 42.577478461,# alias for proton
    H2 = 6.535902854, # ½ 𝑆, 0.01% abund.
    D2 = 6.535902854, # alias for deuteron
    H3 = 45.41483815, # ½ 𝑆, 0% abund (unstable)
    T3 = 45.41483815, # alias for trition
    He3 =-32.43604456,#  ½ 𝑆, 0.001 % abund.
    Li6 = 6.266099,   #  1 𝑆, 7.59% abund.
    Li7 = 16.548177,  # ³/₂𝑆, 92.4% abund.
    C13 = 10.707746,  #  ½ 𝑆, 1.07% abund.
    N14 = 3.076273,   #  1 𝑆, 99.6% abund.
    N15 =-4.315255,   #  ½ 𝑆, 0.36% abund.
    F19 = 40.069244,  #  ½ 𝑆, 100 % abund.
    Si29 =-8.461871,  #  ½ 𝑆, 4.7 % abund.
    P31 = 17.241162,  #  ½ 𝑆, 100 % abund.
    Se77 = 	8.13422,  #  ½ 𝑆, 7.6 % abund.
    Y89 = -2.093134,  #  ½ 𝑆, 100 % abund.
    Rh103 = -1.34600, # ½ 𝑆, 100 % abund.
)
const dpath = joinpath(homedir(), "nmrdata/data/2025")
# Some more concise and more clear type aliases.
const n64 = Int64
const z64 = Complex{n64}
const f64 = Float64
const c64 = Complex{f64}
# abstract and union types using blackbold (set-theoretic) notation
const ℝ = Real
# AbstractVector is compatible with CuArray and Array types of dim k×1.
const ℂ¹ = AbstractVector{c64}
const ℝ¹ = AbstractVector{f64}
const ℕ¹ = AbstractVector{n64}
const ℤ¹ = AbstractVector{z64}
#Type alias for FID data which, in general, is sparse due to ubiquity of zero-filling
const 𝕊¹ = AbstractSparseVector{c64}
# AbstractMatrix is compatible with CuArray and Array types of dim k×j.
const ℂ² = AbstractMatrix{c64}
const ℝ² = AbstractMatrix{f64}
const 𝕊² = AbstractSparseMatrix{c64}
# AbstractArray is compatible with CuArray and Array types of dim k×j×...
const ℂᴺ = AbstractArray{c64}
const ℝᴺ = AbstractArray{f64}
const 𝕊ᴺ = AbstractSparseArray{c64}

"""
    Abstract Supertype of spectrum structs.
"""
abstract type AbstractSpectrum end
# include structs first
include("bruker/Bruker_structs.jl")
#= todo
# include("varian/Varian_structs.jl")
# include("jeol/Jeol_structs.jl")
=#
# FracSpectrum references other spectrum structs as it contains constructors 
# for each of them.
include("FracSpectrum.jl")
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
include("shapedpulses.jl")
try
    if debug
        @info "Debug flag set."
    end
catch 
    global debug = false
end

export BrukerSpectrum, ProcessedSpectrum, dump, # constructors, I/O
       plot, plot!, # plotting
       lsq_analyze, candidates, decompose, # decomposition
       ppmtoindex, hztoindex, ppmtohz, # unit conversion
       zero_fill,
       baseline_correct!, # processing
       limits, chemical_shifts, extract, copy!, # utility functions
       interpolate, resample, # interpolation
       integrate, # integration
       intrng, intrng_data, intrng_indices, intrng_shifts, # integration regions
       debug # debug variable

data = BrukerSpectrum("test/data/1D_test/1")
