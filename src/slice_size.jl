# function to compute the size of a slice from a given gradient amplitude

function get_SliceLength(g, δ; 
                        γ = 42.577478461 , # γₚ₊ / 2π, MHz T-1
                        Gmax = 50, # Gauss cm-1
                        message = true)
    G = g * Gmax * 1e-4 # G cm-1 => T cm-1 
    Γ = γ * G * 1e6 # => Hz cm-1
    l = δ / Γ # => Hz / (Hz cm-1) = cm
    ℓ = round(l, sigdigits = 3) # prettier
    message && println("\n 
            Selects $ℓ cm given relative gradient strength 0 ≤ g ≤ 1 and pulse width in Hz. \n
            Check that γₚ₊ and Gmax $Gmax G cm⁻¹ are appropriate. \n")
    return l;
end 
