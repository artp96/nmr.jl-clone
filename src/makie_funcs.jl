# A makie function to plot a popt array as a heatmap
#
# takes an array of 1-D spectra, a left and right window, 
# the variables from the popt.protocol, and outputs a 
# heatmap in integral or peak intensity mode.
#
using DataFrames
I = Int64

function lineplot_Complex(s :: S; 
procno = false, palette = nothing, kwargs...) where S <: PoptSpectrum
    # get last procno by default
    !procno ? procno = s.default_proc : procno
    fig = Figure()
    ax = Axis(fig[1, 1]; palette = palette, ) 
    lines!(ax, s.procs[procno].re_ft, label = L"\Re")
    lines!(ax, s.procs[procno].im_ft, label = L"\Im")
    return fig
end

"""
    splatted_heatmaps(a <: AbstractArray, dims <: Tuple, vars <: NamedTuple; kwargs ...)
_______________________________________________________________________
kwargs: 
-   colormap = :plasma
-   split_heatmaps = false #!todo: N separate plots
-   save_fig = "" #! todo: save outputs if string nonempty
"""
function splatted_heatmaps(a :: A, dims :: T, vars :: NT; colormap = :plasma) where {
    A <: AbstractArray,
    T <: Tuple,
    NT <: NamedTuple
    }

    colorrange = extrema(a)
    varnames = String.(keys(vars))
    # Need N heatmaps where N is the smallest dimension in vars.
    N = findmin(dims)[2]
    box = reshape(a, dims)
    slices = eachslice(box, dims = N)

    
    z = "$(varnames[N])"
    x, y = filter(xy -> xy != z, varnames)

    n = length(slices)
    if isinteger(sqrt(n)) 
        (grid_w, grid_h) = (Int.(sqrt(n)), Int.(sqrt(n))) 
    else
        (grid_w, grid_h) = (ceil(I, sqrt(n)), floor(I, sqrt(n)))
    end
    fig = Figure(; size = (1600, 1000))
    for col in 1:grid_w, row in 1:grid_h
        idx = (row - 1) * grid_w + col
        idx > n && continue
        
        zvar = round(vars[Symbol(z)][idx], sigdigits = 4)
        ax = Axis(fig[row, col], 
                  title = "$idx: $z = $zvar",
                  aspect = 1)
        x_ax, y_ax = vars[Symbol(x)], vars[Symbol(y)]
        # Makie tries to automatically sort axis labels - fix this
        ax.xreversed = (x_ax[end] < x_ax[1]) 
        ax.yreversed = (y_ax[end] < y_ax[1])
        hm = heatmap!(ax, x_ax, y_ax,
                 slices[idx], 
                 colormap = colormap,
                 colorrange = colorrange
                 )

        # Only keep x-ticks on bottom row
        if row != grid_h
            hidespines!(ax, :b)
            ax.xticksvisible = false
            ax.xlabelvisible = false
        end

        # Only keep y-ticks on left column
        if col != 1
            hidespines!(ax, :l)
            ax.yticksvisible = false
            ax.ylabelvisible = false
        end
        # ⬆ Hide top spine always (optional)
        hidespines!(ax, :t)

        # ➡ Hide right spine always (optional)
        hidespines!(ax, :r)

        cursorpos = Observable(Point2f(0, 0))

        on(events(ax).mouseposition) do pos
            cursorpos[] = pos
        end
        # TODO! add a tooltip on hover
        #tooltip!(hm, cursorpos) do p
        # Convert screen position to data coordinates
        #    xdata, ydata = to_world(ax.scene, p)
        #     xi = findfirst(>(xdata), x_ax)
        #     yi = findfirst(>(ydata), y_ax)
        #     if isnothing(xi) || isnothing(yi)
        #         return ""
        #     end
        #     val = slices[idx][yi, xi]
        #     return "x=$(x[xi]), y=$(y[yi])\nval=$(round(val, digits=3))"
        # end
    end

    Colorbar(fig[:, grid_w + 1], colormap = colormap, colorrange = colorrange, label = "Intensity a.u.")
    # Add global x/y labels using Label()
    Label(fig[grid_h+1, 1:grid_w], "$x", fontsize=18, tellwidth=false, tellheight=true, halign=:center)
    Label(fig[1:grid_h, 0], "$y", fontsize=18, rotation=pi/2, tellwidth=true, tellheight=false, valign=:center)

    resize_to_layout!(fig)
    return fig
    
end

"""
    splatted_heatmaps(a :: A, df :: df; kwargs...) -> figure
    df <: AbstractDataFrame,
    T <: AbstractFloat,
    A <: AbstractArray{T}
_______________________________________________________________________
Function to take a vector or array of values and a dataframe containing 
the Popt array structure and output a series of n heatmaps, where n is 
the smallest dimension of the popt array.

    kwargs: 
    - colormap = :plasma 
"""
function splatted_heatmaps(a :: A, df :: D; kwargs...) where {
    T <: AbstractFloat,
    A <: AbstractArray{T},
    D <: AbstractDataFrame
}

    # Get the array structure from the data frame, and the axis variables
    dims, vars = get_ArrayPoptDims(df)
  
    fig = splatted_heatmaps(a, dims, vars; kwargs...)
    return fig
end

"""
    splatted_heatmaps(df :: D; kwargs...) = splatted_heatmaps(df.Integral, df; kwargs...)

Convenience wrapper to splat the integrals from a popt dataframe.
"""
splatted_heatmaps(df::D; kwargs...) where D<:AbstractDataFrame = splatted_heatmaps(df.Integral, df; kwargs...)

import GLMakie: lines
# some quick wrappers for lines
lines(s::BrukerSpectrum, kind = real, p = s.default_proc, kwargs...) = lines(s[p], kind; kwargs...)

const LINESTYLES = [:solid, :dash,:solid, :dash, :solid, :dash,:solid, :dash]

function lines(p::S, 
               kind::Function = real; 
               colormap = :plasma, 
               title = "", 
               idxs = missing, 
               axis_func = 
               ppmaxis, 
               kwargs... ) where S<:AbstractSpectrum

    @assert in(kind, [real, imag, complex, abs]) "'kind = f::Function' kwarg must refer to one of 'real', 'imag', 'abs' or 'complex'."
    @assert in(axis_func, [ppmaxis, hzaxis, ζaxis, cmaxis]) "'axis_func = f::Function kwarg must refer to one of 'ppmaxis', hzaxis', 'ζaxis' or 'cmaxis'."
    if !ismissing(idxs) 
        @info "Passed indices to lines() plot are expected to be in ppm format."
        idx1, idx2 = ppmtoindex(p, idxs[1]), ppmtoindex(p, idxs[2]) 
    else
        idx1, idx2 = trim_spectrum(p)
    end
    xaxis = axis_func(p)[idx1:idx2] 
    xticks = round.([getindex(xaxis, i) for i ∈ round.(Int, LinRange(1, idx2-idx1, 10) )]; digits = 3)

    fig = Figure()
    ax = Axis(fig[1, 1], title = title, xlabel = string(axis_func)[1:end-4], ylabel = string(kind), xticks = xticks)
    hidespines!(ax)

    # complex mode needs a bit more setup
    if (kind != complex)
        signal = @view kind(p)[idx1:idx2]
    else
        re = @view real(p)[idx1:idx2]
        im = @view imag(p)[idx1:idx2]
        signal = [re im]
    end
     
    for sig in eachcol(signal)
        iter = sig.indices[2]
        lines!(ax, xaxis, sig,
            colormap = colormap,
            color = :blue,
            linestyle = LINESTYLES[iter],
            linewidth = 1

           )
    end
    ax.xticklabelrotation = -π/6
    ax.xticklabelsize = 14
    ax.xticklabelalign = (:left, :center)
    resize_to_layout!(fig)
    return fig
end


# some lines methods for multi-plotting
lines(v::V, kind::Function=real; labels = [s.expno for s in v], kwargs...) where V<:Vector{<:BrukerSpectrum} = lines([s[s.default_proc] for s ∈ v], kind; labels, kwargs...)


function lines(V::P, kind::Function = real;
               colormap = :plasma,
               title = "No Title",
               idxs = missing, 
               axis_func = ppmaxis, 
               labels = collect(1:length(V)), 
               legtitle = "", 
               preproc::Function = identity, 
               kwargs... ) where P<:Vector{<:AbstractSpectrum}

    labels = string.(labels)
    @assert in(kind, [real, imag, complex, abs]) "'kind = f::Function' arg must refer to one of 'real', 'imag', 'abs' or 'complex'."
    @assert in(axis_func, [ppmaxis, hzaxis, ζaxis, cmaxis]) "'axis_func = f::Function kwarg must refer to one of 'ppmaxis', hzaxis', 'ζaxis' or 'cmaxis'."
    
    # ppm axis should be relative to widest spectrum
    max, pno = findmax([size(v) for v ∈ V])

    p = V[pno] 
    idx1, idx2 = trim_spectrum(p)
    if !ismissing(idxs) 
        @info "Passed indices to lines() plot are expected to be a tuple in idxs = (ppm, ppm) format."
        idx1, idx2 = ppmtoindex(p, idxs[1]), ppmtoindex(p, idxs[2]) 
        # no reason not to fix this here.
        idx1, idx2 = idx1 > idx2 ? (idx2, idx1) : (idx1, idx2)
    end
    xaxis = axis_func(p)[idx1:idx2] 
    xticks = round.([getindex(xaxis, i) for i ∈ round.(Int, LinRange(1, idx2-idx1, 10) )]; digits = 3)


    fig = Figure()
    ax = Axis(fig[1, 1], title = title, xlabel = string(axis_func)[1:end-4], ylabel = string(kind), xticks = xticks)
    hidespines!(ax)

    for (i, p) in enumerate(V)

        idx1, idx2 = trim_spectrum(p)
        if !ismissing(idxs) 
            idx1, idx2 = ppmtoindex(p, idxs[1]), ppmtoindex(p, idxs[2]) 
            # no reason not to fix this here.
            idx1, idx2 = idx1 > idx2 ? (idx2, idx1) : (idx1, idx2)
        end

        xaxis = axis_func(p)[idx1:idx2] 

        # complex mode needs a bit more setup
        if (kind != complex)
            preprocessed = preproc(kind(p))
            signal = @view preprocessed[idx1:idx2]
        else
            re = @view real(p)[idx1:idx2]
            im = @view imag(p)[idx1:idx2]
            signal = [re im]
        end
     
        for sig in eachcol(signal)
            
            iter = sig.indices[2]
            lines!(ax, xaxis, sig,
                colormap = colormap,
                linestyle = LINESTYLES[iter],
                linewidth = 1,
                label = labels[i]
               )
        end
    end

    ax.xticklabelrotation = -π/6
    ax.xticklabelsize = 14
    ax.xticklabelalign = (:left, :center)

    fig[1,2] = Legend(fig, ax, legtitle, framevisible = false)


    resize_to_layout!(fig)
    return fig
end

export heatmap_Popt, splatted_heatmaps
