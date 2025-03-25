
using Optimization, ForwardDiff
using OptimizationOptimJL


"""
    ϕ_correct(s::Vector{C}, φ₀ :: F; φ₁ :: F = 0.) where { F <: AbstractFloat, C <: Complex{F}}
    

    Applies a zeroth and optionally a first order phase correction to a Fourier-transformed spectrum.
    Arguments:
    - `re_ft`: Real part of the Fourier-transformed spectrum.
    - `im_ft`: Imaginary part of the Fourier-transformed spectrum.
    - `φ0`: Zero-order phase correction (in degrees).
    - `φ1`: First-order phase correction (in degrees).

    Returns:
    - `Vector{C}`: The corrected real part of the spectrum, ie.
        the endomorphism 𝑓:𝑠 -> 𝑠.
"""
ϕ_correct(r::V, 𝑖::V, ϕ₀ :: T; ϕ₁::T=0.0) where {T <: AbstractFloat, V <: AbstractVector{T}} = ϕ_correct(s = Complex{T}.(r, 𝑖), ϕ₀; ϕ₁)

function ϕ_correct(s::V, ϕ₀; ϕ₁ = 0.) where { C <: Complex, V <: AbstractVector{C} }    
    ϕ₀, ϕ₁ = deg2rad(ϕ₀), deg2rad(ϕ₁);
    n = length(s)
    # normalised frequency axis
    ν = LinRange(0, 1, n) # |> collect |> cu 
    ϕ = exp.(-1im .* (ϕ₀ .+ ϕ₁ .* ν))
    return s .* ϕ
end

# auto phase should take a [[transformed]] spectrum
function auto_ϕ_correct(s :: V; 
            f :: Function = (y -> sum(y[y .< 0.].^2)), 
            ϕ₀:: T = 0., ϕ₁ :: T = 0., tol = false) where {
            T <: AbstractFloat,
            C <: Complex, 
            V <: AbstractVector{C}
            }
    # Ends of the spectrum can distort phase
    idx1, idx2, len = trim_spectra_ends(s)
    s_ = @view s[idx1:idx2]
    # Must meet criteria for OptimizationFunction, AutomaticDifferentiable
    function objective(ϕ, s_)
        ϕ₀, ϕ₁ = ϕ
        y = ϕ_correct(s_, ϕ₀; ϕ₁ = ϕ₁) |> real
        # f can be any loss function, but we seek positive phase for Re(s)
        return f(y)
    end

    ϕ = [ϕ₀, ϕ₁]    
    optF = OptimizationFunction(objective, AutoForwardDiff())
    prob = OptimizationProblem(optF, ϕ, s)
    ϕ_opt = solve(prob, BFGS()) .% 360
    s = ϕ_correct(s, ϕ_opt[1]; ϕ₁ = ϕ_opt[2])
    return s, ϕ_opt
end


import Base.Threads
# function to inplace auto
"""
    auto_ϕ_correct(S <: PoptSpectrum) -> S <: PoptSpectrum 
Phase correct each 1-D fibre of a POPT array serfile, returning a corrected array.
"""
function auto_ϕ_correct(s :: P) where P <: PoptSpectrum#{
#     T <: AbstractFloat,
#     C <: Complex{T},
#     V <: AbstractVector{C},
#     A <: AbstractArray{C}
# }
    fids = eachfibre(s.ser)
    debug && print(eachfibre(fids))
    count = 0
    Threads.@threads for f in fids
        f .= auto_ϕ_correct(f)[1]
        if debug
            lock(ctr_lock)
            count += 1
            println(count)
            unlock(ctr_lock)
        end
    end
    s.ser = collect(fids)
    return s
end

"""
    trim_spectra_ends(S) -> (idx1, idx2, new_length)

Ends of the spectrum are often distorted by DSP, and look ugly in plots or have weird phase errors.
    
"""
function trim_spectra_ends(S; pc = 0.01)
    n = length(S)
    idx1, idx2 = ceil(Int,pc*n), floor(Int,(1-pc)*n)
    return (idx1, idx2, idx1-idx2)
end

export auto_ϕ_correct, ϕ_correct
