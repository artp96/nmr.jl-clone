pkgs = [
 "DataFrames","DSP","ProgressBars","DataFrames","GLMakie","LaTeXStrings","Interpolations","LinearAlgebra","Statistics","RecipesBase","FFTW","FastTransforms","ForwardDiff","GLMakie","Interpolations","LaTeXStrings","LinearAlgebra","Optimization","Optim","ProgressBars","Statistics","StatsBase","RecipesBase","CUDA", "BenchmarkTools"]
for p in pkgs
    Pkg.add(p)
end

using BenchmarkTools
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
