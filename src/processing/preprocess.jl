"""
    preprocess(s<:AbstractSpectrum{T,N}) -> f(t) <: ℂᴺ.
-----------------------------------------------------------------------------
Takes an input spectrum and applies all basic preprocessing to it which under
no conceivable circumstance would be undesirable. That is, 
For Bruker expno (time-domain) data:
    1. The fid is made complex, even values are multiplied by 𝑖.
    2. the GRPDLY circshift is applied in the first dimension, 
    3. any residual ϕ¹ term added to the procpar value
"""
function preprocess(s<:BrukerSpectrum) 
    fid_DSF_undo(s)    


end
