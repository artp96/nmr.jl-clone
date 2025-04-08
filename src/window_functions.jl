# Apply window functions to spectra
# These functions are most simply applied to the FID, then FT
#
#
"""
    exp_window(S, b)
    - S <: AbstractSpectrum
    - b is the broadening factor, in Hz.
-----------------------------------------------------------------------
Apply a gaussian broadening function.
"""


"""
    match_filter(s, f) -> f(s)
    - S <: AbstractSpectrum
    - f(S, b) <: Function 
-----------------------------------------------------------------------
Returns an anonymous function compiled and specialised to exactly the 
spectrum s and matched to the natural linewidth of s which can be simd
or @. efficiently.
"""
