using NMR_acp
using CUDA
using DifferentialEquations, LinearAlgebra, TerminalLoggers
using GLMakie
#CairoMakie.Activate!()
##### THIS HAS BEEN MOVED TO IT'S OWN PROJECT #####

const x = [1., 0., 0.]
const y = [0., 1., 0.]
const z = [0., 0., 1.]
"""
    This file contains a Bloch equation integrator.
    Presently, it works for the simplest system, a single half-spin, 
    and is designed to account for spatial dependency in 𝑧.

"""
abstract type AbstractPulse  end
struct HardPulse <: AbstractPulse
    # timings
    on ::Float64
    off::Float64
    # field during pulse - x, y, z components, ∈ ℝ.
    B :: V where V <: AbstractVector
    # Constructor for type safety
    HardPulse(on, off, B) = new(Float64(on), Float64(off), B)
end

abstract type AbstractInputStruct end
@kwdef mutable struct BlochInputStruct <: AbstractInputStruct
    # Constants
    γ = -2π * 42.58e6      # /s/Tesla for proton
    T₁ = 1.5              # seconds
    T₂ = 0.2              # seconds
    M₀ = [0., 0., 1.0]    # equilibrium magnetization (along z)
    p1= 15e-6             # calibrated 90° pulse duration.
    B₁max = 0.5π/(γ*p1) # Tesla, 
    # Static B field
    B₀ = [0.0, 0.0, 0.0,]  # Default is rotating frame of reference
    # a vector of pulses.
    # by default contains a 15 μs 90° hard pulse. 
    pulprog = [HardPulse(0, p1, B₁max .* x)] 
    # flag variable which is controlled by callbacks outside the differential.
    # when off, B1 is not evaluated. Assume the interesting stuff always starts
    # with a pulse...
    pulse_on::Bool = true
end

@kwdef mutable struct InputStructNoRelaxation <: AbstractInputStruct
    # Constants
    γ = 2π * 42.58e6      # rad/s/Tesla for proton
    M₀ = [0., 0., 1.0]    # equilibrium magnetization (along z)
    p1 = 15e-6            # calibrated 90° pulse duration.
    B₁max = 0.5π /(γ*p1)# Tesla, 
    # Static B field
    B₀ = [0.0, 0.0, 0.0,]  # Tesla, default rotating frame
    # a vector of pulses.
    # by default contains a 15 μs 90° hard pulse. 
    pulprog = [HardPulse(0, p1, B₁max .* x)] 
    # flag variable which is controlled by callbacks outside the differential.
    # when off, B1 is not evaluated. Assume the interesting stuff always starts
    # with a pulse...
    pulse_on::Bool = true
end
default_input = BlochInputStruct()
make_hard_pulse(θ, axis = x, inp=InputStruct()) = HardPulse(0, inp.p1, θ / (inp.γ .* inp.p1) .* axis)
# Magnetic field function (can extend later to include pulses)
pulse_amplitude(pulse::HardPulse, t) = pulse.on ≤ t ≤ pulse.off ? pulse.B : [0., 0., 0.]
# Bloch ODE system
function bloch_rhs!(dM, M, p::BlochInputStruct, t)
    γ, T₁, T₂, M₀ = p.γ, p.T₁, p.T₂, p.M₀
    B = sum([pulse_amplitude(pulse, t) for pulse in p.pulprog]) .+ p.B₀ 
    dM .= γ * cross(M, B)
    dM[1] -= M[1] / T₂
    dM[2] -= M[2] / T₂
    dM[3] -= (M[3] - M₀[3]) / T₁
    # Debugging
    return dM
end
function bloch_rhs!(dM, M, p::InputStructNoRelaxation, t)
    γ, M₀ = p.γ, p.M₀
    B = sum([pulse_amplitude(pulse, t) for pulse in p.pulprog]) .+ p.B₀ 
    dM .= γ * cross(M, B)
    # Debugging
    if 1e-5 < t < 2e-5
    #   @debug "t = $t, dM = $dM, Bfield = $Bfield"
    end
    return dM
end

function plot_trajectory_slider(s::S) where S <: SciMLBase.AbstractTimeseriesSolution
    f = Figure()
    a = Axis3(f[1, 1], title = L"$\mathbf{M}(t)$, rel.")
    ts = s.t
    u₁ = getindex.(s.u, 1)
    u₂ = getindex.(s.u, 2)
    u₃ = getindex.(s.u, 3)
    slider = Slider(f[2, 1], range = eachindex(ts), startvalue = 1)
    lines!(a, u₁, u₂, u₃, #marker = :cross,
            color = ts,
            colormap = :viridis,
            colorrange = extrema(ts),
            transparency = true,
            linewidth = 0.8
             )
    current_idx = slider.value   
    point_obs = @lift(Point3f(u₁[$current_idx], u₂[$current_idx], u₃[$current_idx]))
    scatter!(a, point_obs; markersize = 10, color = :red)
    return f
end

function plot_trajectory_simple(s::S) where S <: SciMLBase.AbstractTimeseriesSolution
    fig = Figure()
    ax = Axis3(fig[1, 1], title = L"$\mathbf{M}(t)$, rel.")
    ts = s.t
    u₁ = getindex.(s.u, 1)
    u₂ = getindex.(s.u, 2)
    u₃ = getindex.(s.u, 3)
    scatter!(ax, u₁, u₂, u₃, #marker = :cross,
            color = ts,
            colormap = [:red, :orange, :brown],
            colorrange = extrema(ts),
            transparency = true,
            markersize = 0.3
            )
    return fig
end

"""
    a function to take a system struct, extract the pulse program timings, set these up as DiscreteCallbacks
    to the solver, force evaluation during these timings, and solve the system
"""
function integrate_bloch_segments(params::AS, tlims = (0., 1.); kwargs...) where AS<:AbstractInputStruct
    pulses = params.pulprog
    savetimes = pulprog_get_timings(pulses)
    #pulse_callbacks = make_pulse_callbacks(pulses) 
    #pulse_edges = unique(vcat(0., [p.on for p in pulses]..., [p.off for p in pulses]...))
    
    prob = ODEProblem(bloch_rhs!, params.M₀, tlims, params)
    sol = solve(
        prob, 
        Tsit5(), 
        abstol = 1e-6, 
        reltol = 1e-6,
        #callback = pulse_callbacks,        
        progress = true,
        tstops = savetimes,
        save_everystep = true, kwargs...
    )
    return sol
end

"""
    A function to parse a list of pulses into a list of timings to force evaluation within the 
    DifferentialEquations solver with saveat=(). Outputs a vector of floats. Optionally specify 
    a precision, default is 1 μs.

"""
function pulprog_get_timings(pulses::V; precision = 1e-7) where {P<:AbstractPulse, V<:Vector{P}}
    return vcat([range(p.on, p.off; step = precision) |> collect for p in pulses]...)
end

"""Fast checking function"""
function in_pulse(t, pulses)
    @inbounds for p in pulses
        if p.on ≤ t ≤ p.off
            return true
        end
    end
    return false
end
#=
"""Function to force B field evaluation of pulses within bloch_rhs!()"""
function make_pulse_callbacks(pulses::V) where {P<:AbstractPulse, V<:Vector{P}}

    # force evaluation of B₁ when t falls within the window of any pulse
    on_cond(u, t, integrator) = in_pulse(t, pulses) && !integrator.p.pulse_on
    function on_affect!(integrator) 
        integrator.p.pulse_on = true
        @debug "pulse ON at t = $(integrator.t)"
    end

    # skip evaluation of B₁ when t is outside the window of any pulse
    off_cond(u, t, integrator) = !in_pulse(t, pulses) && integrator.p.pulse_on
    function off_affect!(integrator)
        integrator.p.pulse_on = false
        @debug "pulse OFF at t = $(integrator.t)"
    end

    return CallbackSet(
        DiscreteCallback(on_cond, on_affect!),
        DiscreteCallback(off_cond, off_affect!)
    )
end
=#

sol = integrate_bloch_segments(default_input)
figslide = plot_trajectory_slider(sol)
figstatic = plot_trajectory_simple(sol)

function plot_trajectory_lines(s::S) where S <: SciMLBase.AbstractTimeseriesSolution
    fig = Figure()
    ax = Axis3(fig[1, 1], title = L"$\mathbf{M}(t)$, rel.")
    ts = s.t
    u₁ = getindex.(s.u, 1)
    u₂ = getindex.(s.u, 2)
    u₃ = getindex.(s.u, 3)
    lines!(ax, u₁, u₂, u₃, #marker = :cross,
            color = ts,
            colormap = [:red, :purple],
            colorrange = extrema(ts),
            transparency = true,
            linewidth = 0.1
            )
    return fig
end
