# module imports
using NMR, GLMakie, CairoMakie

GLMakie.activate!()

import NMR: SFO1, GMAX

const dpath = joinpath(homedir(), "nmrdata/data/2025/b800mib/20250417")
const Dpath = joinpath(homedir(), "nmrdata/data/2025/b800mib/20250423-nug/5")
dosydata = BrukerSpectrum(Dpath)
data = multiimport(dpath); expnos = [d.expno for d in data]
ws = FracSpectrum.(data)
S = (; (Symbol(expnos) .=> data)... )
G = GMAX["CP-TCI-800S4"].z

# ∂Hz / ∂ cm for read and selection gradients
gpz1_Hz = get_ω_dispersion(0.01  ; Gmax = G, message = false);
gpz5_Hz = get_ω_dispersion(0.05  ; Gmax = G, message = false);
gpz10Hz = get_ω_dispersion(0.10  ; Gmax = G, message = false);
gpz20Hz = get_ω_dispersion(0.20  ; Gmax = G, message = false);

# spatial offsets
spoffs2 = LinRange(1.25gpz20Hz, -1.25gpz20Hz, 13); (Hz,cm) = (spoffs2, spoffs2 ./ gpz20Hz)

NOE_whm = get_ActiveVolume(0.2, 7.85; instrument = "MIBB800", Gmax = G, message = false)
horn_Δ  = get_ActiveVolume(0.2, 5.72; instrument = "MIBB800", Gmax = G, message = false)

zgr_vol = get_ActiveVolume(0.2, 7.85; instrument = "MIBB800", Gmax = G, message = false)

# measured bandwidths 
Δδ = [587, 1087, 1928] ./ 8e2; kHz_bandwidth = [5e3, 10e3, 17.2e3]
kHz_bws = NamedTuple( Symbol.(kHz_bandwidth) .=> get_ActiveVolume.(0.01, Δδ; instrument = "MIBB800", Gmax = G, message = false) )
slice_w = get_SliceLength.([0.05, 0.1, 0.2], 17.2e3; Gmax = G);

""" tau_to_p(gₛ, gᵣ, d32, p2, s) -> p ∈ ℕ₀.
    - gₛ : readout gradient intensity
    - gᵣ : selection gradient intensity
    - d32 : gron/groff time
    - p2 : shaped pulse length
    - s : splitting in Hz
--------------------------------------------------------------------
Compute the coherence order p of a given echo refocused by the read
gradient. 
"""
function tau_to_p(gₛ, gᵣ, d32, p2, s)
    # compute gradient areas
    gₛ *= 2d32+p2
    gᵣ *= d32+1/s
    return (gₛ / gᵣ)
end
# calculating this for our 153 Hz modulation:
tau_to_p(20, 1, 300e-6, 51e-6, 153) # ≈ 2.2
#some common plot paramaters
idxs = (-2.659,12.659)
figdir = joinpath(homedir(), "noteforGAM0425", "figures")
# first figure, spectra need reordering
title = "dB Calibration with no selection gradient, 1 % read."
legend = [-6, -7, -8, -9, -10, -11, -11.4]; legtitle = "p2 dB"
plt2data = append!(reverse(data[18:20]), data[21:24])
plot1 = lines(plt2data, labels = legend, idxs = idxs, title = title, legtitle = legtitle)
save(joinpath(figdir, "no_selection_powercal.pdf"), plot1; backend = CairoMakie)


# second figure, with read gradient, small slice, good control?
title = "dB Calibration with 20% selection (3 mm slice), 1 % read."
legend = [-6, -7, -8, -9, -10, -11]; legtitle = "p2 dB"
nums = 27:32
plot1 = lines(data[nums], labels = legend, idxs = idxs, title = title, legtitle = legtitle)
save(joinpath(figdir, "17kHz20pc_powercal.pdf"), plot1; backend = CairoMakie)


