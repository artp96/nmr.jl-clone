# Apply window functions to spectra
# These functions are most simply applied to the FID, then FT
#
#
"""
    Get the fid time indices for each dimension of the fid
"""
function get_times(s::BrukerSpectrum) 
    

end


"""
    apodise!(S, S ∈ ℂᴺ, 𝑓(𝑡₁, 𝑡₂, 𝑡₃...)) -> S ∈ ℂᴺ. 
--------------------------------------------------------------------
    Applies the window function 𝑓 to the spectrum. 
    Window functions are applied to the time domain, before Fourier
    transform, and should be functions of 𝑡 in each dimension.

"""
function apodise!(out::AS, in::AS, f::Function) where AS<:AbstractSpectrum
    # Get the time indices    
    ts = get_times(in)
    


end



export exp_window, match_filter
