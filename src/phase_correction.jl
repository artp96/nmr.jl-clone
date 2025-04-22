
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
ϕ_correct(r::V, 𝑖::V, ϕ₀ :: T, ϕ₁::T=0.0) where {T <: AbstractFloat, V <: AbstractVector{T}} = ϕ_correct(s = Complex{T}.(r, 𝑖), ϕ₀, ϕ₁)

function ϕ_correct(s::V, ϕ₀, ϕ₁) where { C <: Complex, V <: AbstractVector{C} }    
    ϕ₀, ϕ₁ = deg2rad(ϕ₀), deg2rad(ϕ₁);
    n = length(s)
    # normalised frequency axis
    ν = LinRange(0, 1, n) # |> collect |> cu 
    ϕ = exp.(-1im .* (ϕ₀ .+ ϕ₁ .* ν))
    return s .* ϕ
end

"""
    auto_ϕ_correct(s :: V; 
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
-------------------------------------------------------------------------------
Automatic phase correction function to take a 1-D spectrum and apply zeroth and 
first order phase corrections ϕ₀, ϕ₁.
"""
function auto_ϕ_correct(s :: V, idxs; 
            f :: Function = (y -> sum(y[ y .< 0. ] .^ 2)), 
                        ϕ₀:: T = 0.,
                        ϕ₁ :: T = 0.,
                        tol = 1e-3,
                        ) where {
            T <: AbstractFloat,
            C <: Complex, 
            V <: AbstractVector{C}
            }
    idx1, idx2 = idxs
    s_ = @view s[idx1:idx2]
    # Must meet criteria for OptimizationFunction, AutomaticDifferentiable
    function objective(ϕ, s_)
        ϕ₀, ϕ₁ = ϕ
        y = ϕ_correct(s_, ϕ₀, ϕ₁) |> real
        # f can be any loss function, but we seek positive phase for Re(s)
        return f(y)
    end

    ϕ = [ϕ₀, ϕ₁]    
    optF = OptimizationFunction(objective, AutoForwardDiff())
    prob = OptimizationProblem(optF, ϕ, s)
    ϕ_opt = solve(prob, BFGS()) .% 360
    s = ϕ_correct(s, ϕ_opt[1], ϕ_opt[2])
    return s, ϕ_opt
end


import Base.Threads
# function to inplace auto
"""
    auto_ϕ_correct(S0 <: PoptSpectrum) -> S1 <: PoptSpectrum 
_________________________________________________________________________________
Phase correct each 1-D fibre of a POPT array serfile, returning a corrected array.
"""
function auto_ϕ_correct(s :: P, δ::Real; kwargs...) where P <: PoptSpectrum
    idxs = trim_spectrum(s, δ)
    fibres = eachfibre(s.ft * s.ser)
    auto_ϕ_correct(fibres, idxs; kwargs...)
end

"""
    auto_ϕ_correct(fibs :: fibreIterator)
_________________________________________________________________________________
Inner function to iterate over fibres.
"""
function auto_ϕ_correct(fibres :: F, idxs; kwargs...) where {F <: fibreIterator} 
    isempty(idxs) && @error(AssertionError("Empty tuple propagated to fibreIterator auto_ϕ_correct, idxs should be defined at this point."))
    if debug
        println("idxs $idxs")
        count = 0
    end
    
    # multithreading makes a big difference on large arrays.
    Threads.@threads for f in fibres
        
        f .= auto_ϕ_correct(f, idxs; kwargs...)[1]

        if debug
            lock(ctr_lock)
            count += 1
            debug && println(count * "Updated fibre: ", real(f[1:5]))  # Debug view
            unlock(ctr_lock)
        end
    end
    return fibres
end


function auto_ϕ_correct2(fibres :: F, idxs; kwargs...) where {F <: fibreIterator} 
    
    isempty(idxs) && @error(AssertionError(
        "Empty tuple propagated to fibreIterator auto_ϕ_correct, idxs should be defined at this point."
    ))

    if debug
        println("idxs $idxs")
        count = 0
    end
    
    # multithreading makes a big difference on large arrays.
    for f in fibres
        
        f .= auto_ϕ_correct2(f, idxs; kwargs...)[1]


        if debug
            lock(ctr_lock)
            count += 1
            debug && println(count * "Updated fibre: ", real(f[1:5]))  # Debug view
            unlock(ctr_lock)
        end
    end
    return fibres
end

function auto_ϕ_correct2(s :: V, idxs; 
                         L :: Function = _imag_loss, 
                        ϕ₀:: T = 0.,
                        ϕ₁ :: T = 0.,
                        tol = 1e-3,
                        apk0 = true,
                        apk1 = true
                        ) where {
            T <: AbstractFloat,
            C <: Complex, 
            V <: AbstractVector{C}
            }
    idx1, idx2 = idxs
    s_ = @view s[idx1:idx2]

    if apk0    
        f₀(s, ϕ₀) = ϕ_correct(s, ϕ₀, ϕ₁)
        ϕ₀ = _bisection_solver(f₀, s_, L; y0 = ϕ₀, tol = tol)
    end
    
    if apk1
        f₁(s, ϕ₁) = ϕ_correct(s, ϕ₀, ϕ₁)
        ϕ₁ = _bisection_solver(f₁, s_, L; y0 = ϕ₁)
    end

    s = ϕ_correct(s, ϕ₀, ϕ₁)
    return s, (ϕ₀, ϕ₁)
end



import GLMakie: lines
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
function _bisection_solver(f :: F1, X :: A, loss :: F2;
                           tol = 1e-3,
                           y0 = missing,
                           max_iter = 1e3,
                           debug = false) where {
    F1 <: Function, F2 <: Function, A <: AbstractArray{<:Complex} }
    max_iter = Int(max_iter)
    
    # f, l should be def. s.t. f:ℂᴺ -> ℂᴺ, l:ℂᴺ -> ℝ¹
    g = loss ∘ f
    # Initial phase guess - from first period
    # sample the phase correction
    initials = [g(X, ϕ) for ϕ in LinRange(0, 64π, 128)]

    max, idx = findmax(initials)
    εₗ, εᵣ = initials[idx-1], initials[idx+1]
    ϕₗ, ϕᵣ = (idx-1)π, (idx+1)π 

    if debug
        ϕs = zeros(Float64, max_iter)
        εs = zeros(Float64, max_iter)
    end
    # bracket the maximum
    if εᵣ > εₗ
        ϕₗ, εₗ = idx*π, max 
    else
        ϕᵣ, εᵣ = idx*π, max
    end
    iter = 1
    ϕₙ = 0.
    εₙ = 0.
    while iter < max_iter && tol < abs((εᵣ-εₗ)/εᵣ)
        ϕₙ = (ϕᵣ+ϕₗ)/2
        εₙ = g(X, ϕₙ)
        
        # check for maximum in left half
        if εₙ - εₗ < εₙ - εᵣ
            ϕᵣ = ϕₙ
            εᵣ = εₙ
        else
            ϕₗ = ϕₙ
            εₗ = εₙ
        end
        if debug
            ϕs[iter] = ϕₙ
            εs[iter] = εₙ
        end

        iter += 1
    end

    if debug
        filter!(iszero, ϕs)
        filter!(iszero, εs)
    end
    # y_L = ismissing(y0) ? rand() * 2π : (y0 % 2π)
    # debug && println("ln 221: y_L = $y_L")
    # y_R = y_L - π
    #
    # g_L = g(x, y_L) 
    # debug && println("ln 224: g_L = $g_L")
    # g_L isa AbstractFloat ? ε = g_L : throw(error("l ∘ f (x,y) -> ε does not yield a float."))
    #
    # g_R = g(x, y_R)

    
    # error tracking for debugging
    # if debug  
    #     ε_vec = [ε]
    # end
    # mid = y_R - y_L / 2
    # # initial error from left bisector
    #
    # iter = 0
    # while ε > tol && iter < max_iter && abs(y_L - y_R) < 1e-3
    #
    #     mid = y_R - y_L / 2
    #
    #     ε_i = g(x, mid)    
    #     # Seek a plateau; break out of the loop
    #     if ε - ε_i < tol        
    #         break
    #
    #     # handle minima in left half, where ε = ε(y_L)
    #     elseif ε * ε_i < 0.
    #         y_R = mid
    #         # never need to evaluate the error at the right-hand bisector
    #
    #     # handle minima in the right half
    #     else # ε * ε_i > 0.
    #         # cache the error from this mid as the new left hand bisector
    #         y_L = mid 
    #         ε = ε_i
    #     end
    #
    #
    #     # new error is the error of the old mid value
    #     debug && push!(ε_vec, ε_i)
    #     ε = ε_i
    #     iter += 1
    # end
    return debug ? (ϕs, εs) : ϕₙ
end

"""
    _imag_loss(S :: A) where A <: AbstractArray{<: Complex} -> N <: Real
Function to determine the total imaginary component of the fid as a loss function
for phase correction, noting that this should tend to zero as the real part tends to a maximum,
and be pi-periodic in phase.
Complies with l:ℂᴺ -> ℝ¹
"""
_imag_loss(S :: A) where A <: AbstractArray{<: Complex} = -sum(imag(S).^2)
_norm_loss(S :: A) where A <: AbstractArray{<: Complex} = norm(real(S))
L = _imag_loss

export auto_ϕ_correct, ϕ_correct, auto_ϕ_correct2


const test_tols = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
using BenchmarkTools
struct phase_correction_test
    benchmark
    solution
    errors
end
function test_phase_correction(v :: V, idxs :: Tuple, loss) where V <: AbstractVector{<: Complex}# f should be auto_ϕ_correct
    
    idx1, idx2 = idxs
    v_ = @view v[idx1:idx2]
    
    out = Dict()
    pairs = Vector{Any}(undef, length(test_tols))
    for i in eachindex(test_tols)
        tol = test_tols[i]
        f(s, ϕ₀) = ϕ_correct(s, ϕ₀,0.)
        bm = @benchmark _bisection_solver($f, $v_, $loss; tol = $tol) 

        ϕ₀, ε = _bisection_solver(f, v_, loss; tol = tol, debug = true) 


    pairs[i] = phase_correction_test(bm, ϕ₀, ε)

    end
    out = Dict(zip(test_tols,pairs))
    return out     
end
