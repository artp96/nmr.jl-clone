
"""
    Dict of actual Spectrometer ¹H frequencies
"""
const SFO1 = Dict{String, Float64}(
    "B07B500" => 500.1323506,
    "B11B500" => 500.1330883,
    "B10B300" => 299.9914100,
    "B10B400" => 399.9218796,
    "B700" => 700.1343233,
)

"""
    get_SliceLength(g,             # Applied gradient percent, 0 < g ≤ 1,
                    δ;             # Selective window size in Hz;
                    Gmax = 50,     # Maximum gradient strength in G cm⁻¹,
                    message = true # Print a message to the IO.
                    )

function to compute the size of a slice from a given gradient amplitude and selective pulse bandwidth.
""" 
function get_SliceLength(g, δ; 
                        γ = 42.577478461, # γₚ₊ / 2π, MHz T-1
                        Gmax = 50, # Gauss cm-1
                        message = true)
    !(0. ≤ g ≤ 1.) && throw(ArgumentError("Applied gradient g should be a relative value 0 ≤ g ≤ 1.")) 
    if g != 0
        G = g * Gmax * 1e-4 # G cm-1 => T cm-1 
        Γ = γ * G * 1e6 # => MHz cm-1 => Hz cm-1
        l = δ / Γ # => Hz / (Hz cm-1) = cm
        ℓ = round(l, sigdigits = 3) # prettier
        message && println("\n 
            Selects $ℓ cm given relative gradient strength $(100*g)% and pulse width in Hz.
            Check that γₚ₊ and Gmax $Gmax G cm⁻¹ are appropriate.")
        return l;
    else
        l = 0
        message && println("\n 
            Selects $l cm given relative gradient strength $(100*g)% and pulse width in Hz.
            Check that γₚ₊ and Gmax $Gmax G cm⁻¹ are appropriate.")
        return l;
    end
end 

""" 
    get_ActiveVolume(g, Δδ; 
                    γ = 42.577478461, # γₚ₊ / 2π, MHz T-1
                    Gmax = 50, # Gauss cm-1
                    message = true)
    
Compute the active volume from a given spectral dispersion (in ppm) with a given gradient value.
    f(g, Δδ) -> l cm
"""
function get_ActiveVolume(g :: T, Δδ :: T;
                    γ = 42.577478461, # γₚ₊ / 2π, MHz T-1
                    instrument = "B07B500", # default for B07
                    Gmax = 50, # Gauss cm-1
                    message = true) where T <: AbstractFloat
    !(0. ≤ g ≤ 1.) && throw(ArgumentError("Applied gradient g should be a relative value 0 ≤ g ≤ 1.")) 
    B0 = SFO1[instrument] * 1e6 # MHz => Hz
    Δω = 1e-6 * Δδ * B0 # (Hz / Hz) * Hz => Hz
    ∂ₗ = get_ω_dispersion(g; γ = γ, Gmax = Gmax, message = false) # Hz cm-1
    l = Δω / ∂ₗ  # Hz / Hz cm-1 => cm
    ℓ = round(l, sigdigits = 3)
    message && println("\n 
    ∂/∂ₗ = $∂ₗ Hz cm-1
    len = $ℓ cm
    Given relative gradient strength $(100*g)% and peak dispersion $Δδ ppm.

    Check that γₚ₊ and Gmax $Gmax G cm⁻¹ are appropriate.")
    return l;
end

function get_ω_dispersion(g; 
                    γ = 42.577478461, # γₚ₊ / 2π, MHz T-1
                    Gmax = 50, # Gauss cm-1
                    message = true)
    !(0. ≤ g ≤ 1.) && throw(ArgumentError("Applied gradient g should be a relative value 0 ≤ g ≤ 1.")) 

    if g != 0
        G = g * Gmax * 1e-4 # G cm-1 => T cm-1 
        Γ = γ * G * 1e6 # T cm-1 MHz T-1 => MHz cm-1 => Hz cm-1
        𝛤 = round(Γ, sigdigits = 3) # prettier
        message && println("\n 
            ∂/∂ₗ = $𝛤 Hz cm-1 given relative gradient strength $(100*g)% and pulse width in Hz.
            Check that γₚ₊ and Gmax $Gmax G cm⁻¹ are appropriate.")
        return Γ;
    else
        l = 0
        message && println("\n 
            ∂/∂ₓ = $l Hz cm-1 given relative gradient strength $(100*g)% and pulse width in Hz.
            Check that γₚ₊ and Gmax $Gmax G cm⁻¹ are appropriate.")
        return l;
    end
    
end
export get_SliceLength, get_ActiveVolume, get_ω_dispersion
