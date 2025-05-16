using NMR
using GLMakie
using BenchmarkTools
using StatsBase
using Debugger

const dpath = joinpath(homedir(), "nmrdata/data/2025/b300b10/20250430/b300b10")
const imported = multiimport(dpath)
const data = multiwrap(imported)

tsts = Dict(d.src.expno => cpu!(d) for d in data[2:5])

import NMR.full_angle_correction

fig0, p_0, dp0, Dp0 = full_angle_correction(tsts[11]; debug = true)

import NMR.auto_ϕ_correct_3

opt1 = auto_ϕ_correct_3(tsts[11], rand(), rand(), 1e-6)
try1 = lines(ϕ_correct(tsts[11].fs, opt1.minimizer...)|>real)

#fig1 = lines(out)
