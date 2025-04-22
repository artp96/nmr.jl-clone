using NMR
import NMR: SFO1, GMAX
const dpath = joinpath(homedir(), "nmrdata/data/2025/b800mib/20250417")

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

data = multiimport(dpath)
# product of gaussian signal profile with bandwidth w centred at spoffs s₀, 
# a linear gradient signal profile, and a d²g/dz² term
#profile(w, s, s₀ = 0.,) = exp(-( (s-s0)/2.355w )^2 ) * ( (s - s₀) + (s^3 - s₀) )

""" freq_profile(s, g)
    - s <: Spectrum
    - g - applied gradient, default "GPZ9" * GMAX
--------------------------------------------------------------------
"""
function freq_profile(s, g)
    s_ζ = real(s)    
    

    

end


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
tau_to_p(20, 1, 350e-6, 51e-6, 153) # ≈ 2.2

