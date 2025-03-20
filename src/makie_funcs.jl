# A makie function to plot a popt array as a heatmap
#
# takes an array of 1-D spectra, a left and right window, 
# the variables from the popt.protocol, and outputs a 
# heatmap in integral or peak intensity mode.
#
using GLMakie, LaTeXStrings

function lineplot_Complex(s :: S; procno = false, palette = nothing, kwargs...)where S <: Main.NMR.PoptSpectrum{T} where T <: AbstractFloat
    # get last procno by default
    !procno ? procno = s.default_proc : procno
    fig = Figure()
    ax = Axis(fig[1, 1]; palette = palette, ) 
    lines!(ax, s.procs[procno].re_ft, label = L"\Re")
    lines!(ax, s.procs[procno].im_ft, label = L"\Im")
    return fig
end


function splatted_heatmaps(a :: A; vars = missing, colormap = :plasma) where {T <: AbstractFloat, A <: AbstractArray{T}}
    d = size(a); I = Int64
    
    colorrange = extrema(a)
    a = eachslice(a, dims = length(d))
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

export heatmap_Popt, splatted_heatmaps
