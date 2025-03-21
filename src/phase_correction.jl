
using Optim


"""
    ϕ_correct(s::Vector{C}, φ₀ :: F; φ₁ :: F = 0.) where { F <: AbstractFloat, C <: Complex{F}}
    

    Applies a zeroth and optionally a first order phase correction to a Fourier-transformed spectrum.
    Arguments:
    - `re_ft`: Real part of the Fourier-transformed spectrum.
    - `im_ft`: Imaginary part of the Fourier-transformed spectrum.
    - `φ0`: Zero-order phase correction (in degrees).
    - `φ1`: First-order phase correction (in degrees).

    Returns:
    - `Vector{T}`: The corrected real part of the spectrum, ie.
    the endomorphism 𝑓:𝑠 -> 𝑠.
"""
ϕ_correct(r::V, 𝑖::V, ϕ₀ :: T; ϕ₁::T=0.0) where {T <: AbstractFloat, V <: AbstractVector{T}} = ϕ_correct(s = Complex{T}.(r, 𝑖), ϕ₀; ϕ₁)

function ϕ_correct(s::V, ϕ₀ :: F; ϕ₁ :: F = 0.) where { F <: AbstractFloat, C <: Complex, V <: AbstractVector{C} }    
    ϕ₀, ϕ₁ = deg2rad(ϕ₀), deg2rad(ϕ₁);
    n = length(s)
    # normalised frequency axis
    ν = LinRange(0, 1, n) # |> collect |> cu 
    ϕ = exp.(-1im .* (ϕ₀ .+ ϕ₁ .* ν))
    return s .*= ϕ
end

function auto_ϕ_correct(s :: V; 
            f :: Function=(y ->sum(y[y .< 0.].^2)), 
            ϕ₀:: T = 0., ϕ₁ :: T = 0., tol = false) where {
        T <: AbstractFloat, C <: Complex, 
        V <: AbstractVector{C}
            }
    function objective(ϕ :: AbstractVector{W}) where W <: AbstractFloat
        ϕ₀, ϕ₁ = ϕ
        y = ϕ_correct(s, ϕ₀; ϕ₁ = ϕ₁) |> real
        return f(y)
    end

    ϕ = [ϕ₀, ϕ₁]    
    res = optimize(objective, ϕ)
    ϕ_opt = Optim.minimizer(res)
    s = ϕ_correct(s, ϕ_opt[1]; ϕ₁ = ϕ_opt[2])
    return s, ϕ_opt
end

# function to in-place determine the fft from a FID
auto_ϕ_correct(ser :: AbstractVector{Re}) where Re <: Real = auto_ϕ_correct(fft(ser))
# function to inplace auto
function auto_ϕ_correct(s :: S) where S <: PoptSpectrum
    return "WIP"
end

