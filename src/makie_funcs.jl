# A makie function to plot a popt array as a heatmap
#
# takes an array of 1-D spectra, a left and right window, 
# the variables from the popt.protocol, and outputs a 
# heatmap in integral or peak intensity mode.
#

function lineplot_Complex(s :: S; procno = false, palette = nothing, kwargs...) where S <: PoptSpectrum
    # get last procno by default
    !procno ? procno = s.default_proc : procno
    fig = Figure()
    ax = Axis(fig[1, 1]; palette = palette, ) 
    lines!(ax, s.procs[procno].re_ft, label = L"\Re")
    lines!(ax, s.procs[procno].im_ft, label = L"\Im")
    return fig
end

I = Int64

"""
    splatted_heatmaps(a :: A, df :: df; colormap = :plasma) where {
    df <: AbstractDataFrame,
    T <: AbstractFloat,
    A <: AbstractArray{T}
    } -> fig

Function to take a vector or array of values and a dataframe containing the Popt array structure and output a series
of n heatmaps, where n is the smallest dimension of the popt array.
"""
function splatted_heatmaps(a :: A, df :: DataFrame; colormap = :plasma) where {T <: AbstractFloat, A <: AbstractArray{T}}
    d = size(a); 
    colorrange = extrema(a)

    # Get the array structure from the data frame, and the axis variables
    dims, vars = get_ArrayPoptDims(df)

    # Need N heatmaps where N is the smallest dimension in vars.
    N = findmin(dims)[2]
    a = eachslice(a, dims = N)

    n = length(a)
    isinteger(sqrt(n)) ? (w, h) = (Int.(sqrt(n)), Int.(sqrt(n))) : (w, h) = (ceil(I,sqrt(n)), floor(I,sqrt(n)))
    fig = Figure(; size = (1200,800))
    for w in 1:w, h in 1:h
        idx = w * (1 + (h-1))
        ax = Axis(fig[w, h], xlabel = "X (units)", ylabel = "Y (units)", title = "$idx (units)")
        heatmap!(ax, a[w + h - 1], colormap=colormap, colorrange = colorrange)
    end
    Colorbar(fig[:, h+1], colormap = colormap, colorrange = colorrange, label = "Intensity a.u.")
    resize_to_layout!(fig)
    return fig
end

"""
    splatted_heatmaps(df :: D; kwargs...) = splatted_heatmaps(df.Integral, df; kwargs...)

Convenience wrapper to splat the integrals from a popt dataframe.
"""
splatted_heatmaps(df :: D; kwargs...) where D <: AbstractDataFrame = splatted_heatmaps(df.Integral, df; kwargs...)

export heatmap_Popt, splatted_heatmaps
