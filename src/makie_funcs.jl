# A makie function to plot a popt array as a heatmap
#
# takes an array of 1-D spectra, a left and right window, 
# the variables from the popt.protocol, and outputs a 
# heatmap in integral or peak intensity mode.
#
using GLMakie, LaTeXStrings

function lineplot_Complex(s :: S; procno = false, palette = nothing, kwargs...) where S <: Main.NMR.Spectrum
    # get last procno by default
    !procno ? procno = s.default_proc : procno
    fig = Figure()
    ax = Axis(fig[1, 1]; palette = palette, ) 
    lines!(ax, s.procs[procno].re_ft, label = L"\Re")
    lines!(ax, s.procs[procno].im_ft, label = L"\Im")
    return fig
end
