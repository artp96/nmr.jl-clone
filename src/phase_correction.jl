
using Optimization, ForwardDiff
using Base: debug_color
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
                        ϕ₀:: T = 0.,
                        ϕ₁ :: T = 0.,
                        tol = 1e-3,
                        pivot_ppm = false
                        ) where {
            T <: AbstractFloat,
            C <: Complex, 
            V <: AbstractVector{C}
            }
    idx1, idx2, len = trim_spectrum(s; pivot = pivot_ppm)
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
    auto_ϕ_correct(S0 <: PoptSpectrum) -> S1 <: PoptSpectrum 
Phase correct each 1-D fibre of a POPT array serfile, returning a corrected array.
"""
function auto_ϕ_correct(s :: P; kwargs...) where P <: PoptSpectrum#{
#     T <: AbstractFloat,
#     C <: Complex{T},
#     V <: AbstractVector{C},
#     A <: AbstractArray{C}
# }
    fids = auto_ϕ_correct(eachfibre(s.ser); kwargs...)
    #! TODO there's a type error here 260325
    s.ser = collect(fids)
    return s
end

"""
    auto_ϕ_correct(fibs :: fibreIterator)
Inner function to iterate over fibres.
"""
function auto_ϕ_correct(fibres :: F; kwargs...) where {F <: fibreIterator} 
    
    if debug
        count = 0
    end
    # multithreading makes a big difference on large arrays.
    Threads.@threads for f in fibres
        
        fft!(f)
        f .= auto_ϕ_correct(f; kwargs...)[1]
        if debug
            lock(ctr_lock)
            count += 1
            debug && println(count)
            unlock(ctr_lock)
        end
        ifft!(f)
    end
    return fibres
end

"""
    trim_spectrum(S) -> (idx1, idx2, new_length)

Ends of the spectrum are often distorted by DSP, and look ugly in plots or have weird phase errors.
    
"""
function trim_spectrum(S; pc = 0.01, pivot = false)
    # default trim
    if !pivot
        n = length(S)
        idx1, idx2 = ceil(Int,pc*n), floor(Int,(1-pc)*n)
        return (idx1, idx2, idx1-idx2)
    end
end




"""
    _bisection_solver(f :: F1, x :: A, tol, loss :: F2; y0 = missing, max_iter = 1e3) where {
    F1 <: Function, F2 <: Function, A <: AbstractArray{<:Complex} } -> y_opt
    
Implements the bisection method rootfinding algorithm seeking the zero of a loss function 
composed with a function of some data x ∈ ℂᴺ and a minimisation variable y ∈ ℝ¹. Requires
    f :: F1 f(x, y) : x -> x'∈ ℂᴺ,
    l :: F2 loss(x) : x'-> ε ∈ ℝ¹,
and returns the minimising y value. The composition l ∘ f must be defined to bracket zero 
at π intervals.

Warning: will miss optima if there are an odd number of roots on y ∈ {0, π}, and fail if there are an even number!
"""
function _bisection_solver(f :: F1, x :: A, tol, loss :: F2;
                           y0 = missing, max_iter = 1e3) where {
    F1 <: Function, F2 <: Function, A <: AbstractArray{<:Complex} }
    max_iter = Int(max_iter)
    
    # f, l should be def. s.t. f:ℂᴺ -> ℂᴺ, l:ℂᴺ -> ℝ¹
    g = loss ∘ f
    # Initial phase guess
    y_L = ismissing(y0) ? rand() * 2π : (y0 % 2π)
    y_R = y_L - π

    g_L = g(x, y_L) 
    g_L isa AbstractFloat ? ε = g_L : throw(error("l ∘ f (x,y) -> ε does not yield a float."))

    g_R = g(x, y_R)

    if g_R * g_L > 0
        error("The function must have opposite signs at π intervals. \n
              Check for an even number of roots on y ∈ {0, π}, including zero.")
    end
    
    # error tracking for debugging
    if debug  
        ε_vec = Float64[]
    end
    mid = y_R - y_L / 2
    # initial error from left bisector
    
    iter = 0
    while ε > tol && iter < max_iter && abs(y_L - y_R) < 1e-3
        
        mid = y_R - y_L / 2
   
        ε_i = g(x, mid)    
        # Seek a plateau; break out of the loop
        if ε - ε_i < tol        
            break

        # handle minima in left half, where ε = ε(y_L)
        elseif ε * ε_i < 0.
            y_R = mid
            # never need to evaluate the error at the right-hand bisector

        # handle minima in the right half
        else # ε * ε_i > 0.
            # cache the error from this mid as the new left hand bisector
            y_L = mid 
            ε = ε_i
        end

        
        # new error is the error of the old mid value
        debug && push!(ε_vec, ε_i)
        ε = ε_i
        iter += 1
    end
    return debug ? (mid, ε_vec) : mid
end

"""
    _imag_loss(S :: A) where A <: AbstractArray{<: Complex} -> N <: Real
Function to determine the total imaginary component of the fid as a loss function
for phase correction, noting that this should tend to zero as the real part tends to a maximum,
and be pi-periodic in phase.
Complies with l:ℂᴺ -> ℝ¹
"""
_imag_loss(S :: A) where A <: AbstractArray{<: Complex} = sum(imag(S).^2)
_norm_loss(S :: A) where A <: AbstractArray{<: Complex} = norm((imag(S)))
L = _imag_loss

export auto_ϕ_correct, ϕ_correct
