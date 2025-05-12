using NMR, GLMakie
using BenchmarkTools
g, Gmax = 0.2, 67.
fpath  = joinpath(homedir(), "nmrdata/data/2025/b500b07/250411")

S = [Spectrum(joinpath(fpath, s), 1) for s in string.(10:20)]
s = S[1]
p = s[1]

# zgrg spectrum

# get the z axis in Hz and cm
hz_axis = hzaxis(s)
ζ_axis = hztoζ(hz_axis, g*Gmax)

bandwidth_kHz = [(1e-6*s["P"][3])^-1 for s ∈ S] * 1e-3
bandwidth_cm = [get_SliceLength(g, 1e3w; Gmax = Gmax, message = false) for w in bandwidth_kHz]
bw_labels = string.(round.(bandwidth_kHz; sigdigits = 2) .* " kHz, " .* string.(round.(bandwidth_cm;sigdigits=2)) .* " cm."
