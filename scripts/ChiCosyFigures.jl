include("../src/package.jl")

path08 = joinpath(dpath, "b700/20251208")
path10 = joinpath(dpath, "b700/20251210")

chicosyqf = BrukerSpectrum(path08 * "/5")
cosygpppqf = BrukerSpectrum(path08 * "/6")
# non-paramagnetic sample spectra, 1D and a Jres.
F19_1dzg = BrukerSpectrum(path10 * "/6")
Chorus1d = BrukerSpectrum(path10 * "/7")
ChorusJR = BrukerSpectrum(path10 * "/8")

# paramagnetic spectra
#ChiCo = BrukerSpectrum(path10 * "/14")


