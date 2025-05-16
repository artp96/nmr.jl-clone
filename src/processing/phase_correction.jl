import Base: debug_color
import GLMakie: lines, lines!
export auto_ϕ_correct, ϕ_correct, ϕ_correct!, auto_ϕ_correct2
#! TODO: Automatic phase correction based on entropy minimisation, 
# Reference:
# Chen, Weng, Goh & Garland, 2002: viz. https://doi.org/10.1016/S1090-7807(02)00069-1
# Solve with the simplex method, noting
# 0 ≤ ϕ₀ ≤ 360°
# 0 ≤ ϕ₁ / point  ≤ 360°

"""
    ϕ_correct(s::Vector{C}, ϕ₀;:F=0., ϕ₁::F = 0.) where {F<:AbstractFloat,C<: Complex{F}}
------------------------------------------------------------------------------------------- 
Applies a zeroth and a first order phase correction to a Fourier-transformed spectrum.
Arguments:
- `re_ft`: Real part of the Fourier-transformed spectrum.
- `im_ft`: Imaginary part of the Fourier-transformed spectrum.
OR 
- `s` : Complex spectrum
AND
- `φ0`: Zero-order phase correction (in degrees).
- `φ1`: First-order phase correction (in degrees).
Returns:
- `Vector{C}`: ie. the endomorphism 𝑓:𝑠 -> 𝑠.
"""
function ϕ_correct(s::V, ϕ₀=0., ϕ₁=0.) where { C <: Complex, V <: AbstractVector{C} }    
    return s .* exp.(-1im .* (deg2rad(ϕ₀) .+ (deg2rad(ϕ₁) .* eachindex(s)) ./ length(s)) )
end
function ϕ_correct!(s::V, ϕ₀=0., ϕ₁=0.) where { C <: Complex, V <: AbstractVector{C} } 
    return s .*= exp.(-1im .* (deg2rad(ϕ₀) .+ (deg2rad(ϕ₁) .* eachindex(s) ./ length(s))) )
end
ϕ_correct(r::V, 𝑖::V, ϕ₀ :: T, ϕ₁::T=0.0) where {T <: AbstractFloat, V <: AbstractVector{T}} = ϕ_correct((r.-im*𝑖), ϕ₀, ϕ₁)
#=
#Native GP version-seems to be of no benefit, needless complication, 
#eachindex runs quickly on CuArray types 5/5/25
import CUDA: @cuda, happy
function ϕ_kernel(v, ϕ₀, ϕ₁, len)
    y = threadIdx().y 
    x = threadIdx().x
    idx = x + x * (y - 1)
    v[idx] = ϕ₀ + (ϕ₁ * idx / len)
    return nothing # void
end

function ϕ_correct(s::V, ϕ₀=0., ϕ₁=0.) where { C <: Complex, V <: CuVector{C} }    
    len = length(s)
    ϕ₀ = deg2rad(ϕ₀)# |> CuArray
    ϕ₁ = deg2rad(ϕ₁)# |> CuArray
    ϕs = similar(s)
    @cuda threads = 1024 ϕ_kernel(ϕs, ϕ₀, ϕ₁, len)
    return s .* exp.(-1im .* ϕs)
end
=#
"""
    auto_ϕ_correct(s :: V, idxs = (idx1, idx2); kwargs...) -> s; ϕ₀, ϕ₁; optimization
-------------------------------------------------------------------------------------------
Automatic phase correction function to take a 1-D spectrum and apply zeroth and 
first order phase corrections ϕ₀, ϕ₁. 
Kwargs:
- cost=Function : Loss function for optimisation, must be defined such that
f {|v ∈ ℂᴺ, p ∈ ℜ¹, f(v,p)| ->|ε ≥ 0 ∈ ℝ¹|}, where p is a regularisation parameter (if used) 
- ϕ₀::T <: AbstractFloat :
- ϕ₁::T <: AbstractFloat :
 p=1e-5 <: AbstractFloat : regularisation parameter, C must accept p but may ignore it.
 force_global = false    :  
"""
function auto_ϕ_correct(s :: V, idxs; 
            cost :: Function = _imag_loss,
                         p = 1e-5,
                         ϕ₀:: T = rand(),
                         ϕ₁ :: T = rand(),
                        ) where {
            T <: AbstractFloat,
            C <: Complex, 
            V <: AbstractVector{C}
            }
    idx1, idx2 = idxs
    ψ = @view s[idx1:idx2]
    init = [ϕ₀, ϕ₁] 
    f(init) = -sum(cost(ψ, init[1], init[2], p))
    return optimize(f, init, NelderMead(); kwargs...)
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
function auto_ϕ_correct0(fibres :: F, idxs; kwargs...) where {F <: fibreIterator} 
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
    f :: F1 := f(x, y) : x -> x'∈ ℂᴺ,
    l :: F2 := loss(x) : x'-> ε ∈ ℝ¹,
and returns the minimising y value. The composition l ∘ f must be defined to bracket zero 
at π intervals.
Warning: will miss optima if there are an odd number of roots on y ∈ {0, π}, and fail if there are an even number!
"""
function _bisection_solver(f :: F1, X :: A, loss :: F2; tol = 1e-3, y0 = missing, max_iter = 1000) where {
    F1 <: Function, F2 <: Function, A <: AbstractArray{<:Complex} }
    max_iter = Int(max_iter)
    
    # f, l should be def. s.t. f:ℂᴺ -> ℂᴺ, l:ℂᴺ -> ℝ¹
    g = loss ∘ f
    # Initial phase guess - from first period
    # sample the phase correction
    initials = [g(X, ϕ) for ϕ in LinRange(0, 2π, 10)]

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
    return debug ? (ϕs, εs) : ϕₙ
end

"""
    _imag_loss(S :: A) where A <: AbstractArray{<: Complex} -> N <: Real
-----------------------------------------------------------------------------------------------
Function to determine the total imaginary component of the fid as a loss function
for phase correction, noting that this should tend to zero as the real part tends to a maximum,
and be pi-periodic in phase.
Complies with l:ℂᴺ -> ℝ¹
"""
_imag_loss(S :: A) where A <: AbstractArray{<: Complex} = -sum(imag(S).^2)
_norm_loss(S :: A) where A <: AbstractArray{<: Complex} = norm(real(S))
L = _imag_loss



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
        f(s, ϕ₀) = ϕ_correct(s, ϕ₀, 0.)
        bm = @benchmark _bisection_solver($f, $v_, $loss; tol = $tol) 

        ϕ₀, ε = _bisection_solver(f, v_, loss; tol = tol) 


    pairs[i] = phase_correction_test(bm, ϕ₀, ε)

    end
    out = Dict(zip(test_tols,pairs))
    return out     
end

""" nonneg(A <: AbstractArray) -> f64
-----------------------------------------------------------------------------------------------
Non-negative penalty function. Sum of squared non-negative values, Cu-valid implementation.
(Also semi-optimal for Array).
Default penalty function for spectral_entropy() maximiser.
"""
nonneg(A) = sum( ifelse.(A .< 0, .^(A, 2), 0.) ) / sum(A.^2)

""" spectral_entropy(ψ::A, ϕ₀, ϕ₁, p = 1e-6, P = nonneg) -> f64 
-----------------------------------------------------------------------------------------------
Compute the normalised entropy of the real part of an NMR spectrum given a pair of phase 
offsets.

"""
function spectral_entropy(ψ::A, ϕ₀, ϕ₁, p = 1e-6, P = nonneg) where 
    {
        ℂ <: Complex, 
        A <: AbstractArray{ℂ}
    }
    k = size(ψ, 1)
    ψ = gpu(ψ)
    # setup initial ϕ guesses
    # eq. 6-ii in ref. Chen, Weng, Goh & Garland, JMR. 2002
    R = real( ϕ_correct(ψ, ϕ₀, ϕ₁) )
    hs = zeros(eltype(R), (size(ψ)..., 4) ) |> gpu
    # Derivative via five-point stencil method, assuming periodicity around the edges (usually true)
    slices = enumerate(eachslice(hs; dims = ndims(hs)))
    for (m,h) in slices
        m == 1 ? h .= derivative(R) : h .= derivative(slices.itr[m-1]) 
    end
    # normalise each spectrum
    for h in slices.itr
        h .= abs.(h)
        h ./= sum(h)
    end
    # return entropy
    E = - sum( hs .* log.(hs)) 
    p = p * P(R) / sum(R)
    return E, p
end

    
"""     derivative(V, u = 1)
---------------------------------------------------------------------------------------------
Compute the derivative of V at each point using the five-point stencil method.
Wraps around the boundary, ie. assumes periodicity. Trim the boundary in this case.
The variable u corresponds 
Conserves dimensionality, currently:
    𝑣 ∈ 𝑉 -> ∂ᵤ𝑣 ∈ 𝑉

"""
derivative(v::V, u = 1) where V <: AbstractVector = begin
    ∂v = (-circshift(v, +2) .+ 8 .* circshift(v, +1) .- 8 .* -circshift(v, -1) .- circshift(v, +2)) ./ 12u
    # zero NaN and Inf values arising from zeros in the denominator
    # there is probably a better way to do this but this works for the contexts of 
    # summing the derivative 
    ∂v = ifelse.(-Inf.<∂v.<Inf, ∂v, 0.)
    return ∂v
end 
#! TODO: add a 2-, and 3- or n-D stencil methods


function auto_ϕ_correct_3(ψ::A, ϕ₀, ϕ₁, p, cost = spectral_entropy; kwargs...) where
    {
        ℂ <: Complex,
        A <: AbstractVector{ℂ}
    }
    init = [ϕ₀, ϕ₁] 
    f(init) = -sum(cost(ψ, init[1], init[2], p))
    return optimize(f, init, NelderMead(); kwargs...)
end


function test_correction(d, p = [1e-6, 1e-7, 1e-8, 1e-9]; kwargs...)
    out_nosum  = [auto_ϕ_correct_3(
        fftshift(fft(circshift(d.ft,-d["GRPDLY"]))), 360rand(), 360rand(), p; kwargs...)
        for p in p]

    fig=Figure();
    ax=Axis(fig[1, 1]);
    for k in eachindex(out_nosum)
        spec = ϕ_correct( fftshift( fft(circshift(d.ft, -d["GRPDLY"])) ),
                out_nosum[k].minimizer[1], 
                out_nosum[k].minimizer[2] ) |> cpu
        lines!(ax, real(spec), label = string(p[k]), linewidth = 1) 
    end
    fig[1,2] = Legend(fig, ax, "p weight", framevisible = false)
    resize_to_layout!(fig)
    return fig, out_nosum
end
using GPUArraysCore:@allowscalar
function angle_zero_correct(ψ::A) where {C<:Complex,A<:AbstractArray{C}}
    max, idx = findmax(abs.(ψ))
    @allowscalar ang = angle(ψ[idx])
    return ψ .* exp(-im*(pi/2 - ang))
end

import DSP: unwrap, unwrap!
rel(v) = cpu(v)[10:end-10]./maximum(abs.(v)) |> real |> cpu
# convenience method
full_angle_correction(w::W; kwargs...) where W<:WrappedSpectrum = full_angle_correction(w.fs, w["SW"], w["SFO1"], w["O1"]; kwargs...)

"""A function to calculate PHC0 in ° and PHC1, in ° Hz⁻¹. """
function full_angle_correction(ψ::A, sw, sf, o1; debug = false) where {C<:Complex,A<:AbstractArray{C}}
    
    # get the spectral width in Hz, ppm⋅MHz / len, -> Hz idx⁻¹.
    sw *= 2sf / length(ψ)
    # get the o1p frequency index
    o1p = o1 / sw
    
    # ϕ = ϕ₀ + (ν-ν₀)⋅ϕ₁, solve for ϕ₀, ϕ₁
    # !TODO USE UNWRAP
    # ϕ₀ is global, it arises due to an offset between ϕᵣ and the signal phase.
    # We shall calculate a ϕ₁ first, then ϕ₀ using the residual.
    ϕ = rad2deg.(unwrap(angle.(ψ) .- π |> cpu))
    # get weights as the real signal in the data, 
    w = abs.(ψ)
    # find the pivot ν₀ as the point where the cumulative phase offset is greatest
    # - after this point, the first order correction changes sign
    phmax, ν₀ = findmax(abs.(ϕ))
    
    # get the max with sign
    ϕ₁ = diff(ϕ)
    Φ₁ = derivative(ϕ)

    # a debugging section to plot the phase system 
    if debug
        # make the plot readable
        f=Figure(); ax = Axis(f[1, 1], title = "Phase angle distribution plot (rel.)")
        lines!(ax, rel(real(ψ)), label = "ψ(s)", linewidth = 1)
        lines!(ax, rel(ϕ), label = "|ϕ|", linestyle = :dash, linewidth = 1)
        lines!(ax, rel(cpu(ϕ₁) ./ cpu(ϕ)[2:end]), label = "Δϕ / |ϕ|", linestyle = :dash, linewidth = 1)
        lines!(ax, rel(ϕ₁), label = "Δϕ", linestyle = :dash, linewidth = 1)
        lines!(ax, filter!( ϕₖ -> -Inf < ϕₖ < Inf, rel(Φ₁ |> cpu) ), label = "∂ϕ", linestyle = :dash, linewidth = 1)
        f[1, 2] = Legend(f, ax, "feature", framevisible = false)
        resize_to_layout!(f)
        return f, ϕ, ϕ₁, Φ₁
    end
    return ϕ, ϕ₁
end


