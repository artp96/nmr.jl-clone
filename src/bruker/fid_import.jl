""" fid_DSF_undo(s <: BrukerSpectrum, f = Id) -> fid
-------------------------------------------------------------------------------------------
A special function to apply the appropriate compensation to "undo" an aspect of 
Bruker's DSP which appends data prior to the zero-time point of acquisition for the 
purposes of convolution padding.
Must be called _after_ the FID is made complex. 
Optionally, pass a function f as the second argument to be performed on the FID _without_ 
the prepended DSP convolutional padding. This is valid for e.g. line_broadening(s, hz) or 
zero_fill(s, k) in the Bruker context.

References
 1. W. M. Wrestler and F. Abildgaard, 1996
 2. M. Nilsson & the UoM GNAT collaboration, 2025
"""
function fid_DSF_undo(s :: BrukerSpectrum)
    @assert s.α == 0 "Spectrum does not appear to be time-domain, α = $(s.α)!"
    k = floor(Int, get_DSF_offset(s))
    if iszero(k+ϕ) 
        @warn "Expno $(s.expno) GRPDLY evaluated to zero, 
        \nassuming no DSP to correct!"
        return s
    end
    v = s.fid
    #for fid in eachfibre(fids)
    dsp = deepcopy(v[1:k])
        v = v[k+1:end]
        append!(dsp, v)
    #end
    return v
end

"""
    Function to reapply the effect of Bruker's DSP filter to a FID.
    necessary for applying some time-domain filters.
"""
function fid_DSF_redo(s :: BrukerSpectrum)
    @assert s.α == 0 "Spectrum does not appear to be time-domain, α = $(s.α)!"
    k = floor(Int, get_DSF_offset(s))
    if iszero(k+ϕ) 
        @warn "Expno $(s.expno) GRPDLY evaluated to zero, 
        \nassuming no DSP to correct!"
        return s
    end
    dims = ndims(s.fid)
    v = copy(fid)
    circshift!(s.fid, v, (k,1))
end

""" getoffset(s <: BrukerSpectrum) -> x ∈ u64, y ∈ f64
-------------------------------------------------------------------------------------------
Get the GRPDLY offset or look it up, and return the offset and first-order ϕ₁(ω)
factor in (° Hz⁻¹).
"""
function get_DSF_offset(s :: BrukerSpectrum)::n64
    # simplest case - return GRPDLY iff ∃d > 0.
    if haskey(s, "GRPDLY")
        offset = 0 < s["GRPDLY"] ? s["GRPDLY"] : Throw(ErrorException(
        "GRPDLY defined ≤ 0!
        \nCheck that the values for expno $(s.expno) are listed in the table in nmr.jl/src/fid_import.jl.
        \nGRPDLY = $(haskey(s, GRPDLY) ? s["GRPDLY"] : "missing")
        \nDECIM  = $(haskey(s, DECIM) ? s["DECIM"] : "missing")
        \nDSPFVS  = $(haskey(s, DSPFVS) ? s["DSPFVS"] : "missing")"))
    else
        # otherwise
        row = findfirst(isequal(s["DECIM"]), DSP_TABLE)
        col = s["DSPFVS"] - 9 # Only the 10, 11, 12 cases are listed, why be fancy?
        try 
            offset = DSP_TABLE[row,col]
        catch err
            Throw(ErrorException("Couldn't lookup DSP value!
            \nCheck that the values for expno $(s.expno) are listed in the table in src/fid_import.jl.
            \nGRPDLY = $(haskey(s, GRPDLY) ? s["GRPDLY"] : "missing")
            \nDECIM  = $(haskey(s, DECIM) ? s["DECIM"] : "missing")
            \nDSPFVS  = $(haskey(s, DSPFVS) ? s["DSPFVS"] : "missing")
            \n$err"))
        end
    end

    return -offset
end

const DSP_TABLE=[
#   DECIM      DSPVS10      DSPVS11      DSPVS12
        2       44.750       46.000       46.311;
        3       33.500       36.500       36.530;
        4       66.625       48.000       47.870;
        6       59.083       50.167       50.229;
        8       68.563       53.250       53.289;
       12       60.375       69.500       69.551;
       16       69.531       72.250       71.600;
       24       61.021       70.167       70.184;
       32       70.016       72.750       72.138;
       48       61.344       70.500       70.528;
       64       70.258       73.000       72.348;
       96       61.505       70.667       70.700;
      128       70.379       72.500       72.524;
      192       61.586       71.333           0.;
      256       70.439       72.250           0.;
      384       61.626       71.667           0.;
      512       70.470       72.125           0.;
      768       61.647       71.833           0.;
      1024      70.485       72.063           0.;
      1536      61.657       71.917           0.;
      2048      70.492       72.031           0.]

"""
split_fid(fid :: V{ℝ}) where {V <: AbstractVector, ℝ<:Real} 
Split a fresh FID into it's real and complex parts.
N-D spectra are still acquired time-domain sequentially, so this should work for all N-D spectra.
This function effectively chops the FID length in half, zero filling must be applied 
immediately to get back the original dims.
"""
function split_fid(fid :: V) where V <: ℝ¹
    fid = complex(fid)
    v_re = @view fid[1:2:end]
    v_im = @view fid[2:2:end]
    v_re .-= 1im * v_im
    return v_re
end

export fid_DSF_undo, get_DSF_offset, split_fid


#=
""" APPENDIX
    Based on the procedure described below - A.C.P. 2025.
    Sourced via a comment in UoM GNAT, M Nilsson.
-----------------------------------------------------------------------------------
W. M. Westler and F.  Abildgaard
    July 16, 1996

    The introduction of digital signal processing by Bruker in their DMX
    consoles also introduced an unusual feature associated with the data. The
    stored FID no longer starts at its maximum followed by a decay, but is
    prepended with an increasing signal that starts from zero at the
    first data point and rises to a maximum after several tens of data points.
    On transferring this data to a non-Bruker processing program such as FELIX,
    which is used at NMRFAM, the Fourier transform leads to an unusable spectrum
    filled with wiggles. Processing the same data with Bruker's Uxnmr
    program yields a correct spectrum. Bruker has been rather reluctant
    to describe what tricks are implemented during their processing protocol.

    They suggest the data should be first transformed in Uxnmr and then inverse
    transformed, along with a GENFID command, before shipping the data to another
    program. Bruker now supplies a piece of software to convert the digitally
    filtered data to the equivalent analog form.
    We find it unfortunate that the vendor has decided to complicate
    the simple task of Fourier transformation. We find that the procedure
    suggested by Bruker is cumbersome, and more so, we are irritated since
    we are forced to use data that has been treated with an unknown procedure.
    Since we do not know any details of Bruker's digital filtration procedure
    or the "magic" conversion routine that is used in Uxnmr, we have been forced
    into observation and speculation. We have found a very simple, empirical
    procedure that leads to spectra processed in FELIX that are identical,
    within the noise limits, to spectra processed with Uxnmr. We deposit
    this information here in the hope that it can be of some
    use to the general community.
    The application of a nonrecursive (or recursive) digital filter to time
    domain data is accomplished by performing a weighted running average of
    nearby data points. A problem is encountered at the beginning of
    the data where, due to causality, there are no prior values. The
    weighted average of the first few points, therefore, must include data
    from "negative" time. One naive procedure, probably appropriate to NMR
    data, is to supply values for negative time points is to pad the data with
    zeros. Adding zeros (or any other data values) to the beginning of
    the FID, however, shifts the beginning of the time domain data (FID) to
    a later positive time. It is well known that a shift in the time
    domain data is equivalent to the application of a frequency-dependent,
    linear phase shift in the frequency domain. The 1st order phase shift
    corresponding to a time shift of a single complex dwell point is 360 degrees
    across the spectral width. The typical number of prepended points
    found in DMX digitally filtered data is about 60 data points (see below),
    the corresponding 1st order phase correction is ~21,000 degrees.
    This large linear phase correction can be applied to the transformed data
    to obtain a normal spectrum. Another, equivalent approach is to time
    shift the data back to its original position. This results in the need
    of only a small linear phase shift on the transformed data.
    There is a question as what to do with the data preceding the actual
    FID. The prepended data can be simply eliminated with the addition
    of an equal number of zeros at the end of the FID (left shift). This
    procedure, however, introduces "frowns" (some have a preference to refer
    to these as "smiles") at the edge of the spectrum. If the sweep
    width is fairly wide this does not generally cause a problem. The
    (proper) alternative is to retain this data by applying a circular left
    shift of the data, moving the first 60 or so points (see recommendations
    below) to the end of the FID. This is identical to a Fourier transformation
    followed by the large linear phase correction mentioned above. The
    resulting FID is periodic with the last of the data rising to meet the
    first data point (in the next period). Fourier transform of this
    data results in an approximately phased spectrum. Further linear
    phase corrections of up to 180 degrees are necessary. A zero fill applied
    after a circular shift of the data will cause a discontinuity and thus
    introduce sinc wiggles on the peaks. The usual correction for DC
    offset and apodization of the data, if not done correctly, also results
    in the frowns at the edges of the spectrum.

    In our previous document on Bruker digital filters, we presented deduced
    rules for calculating the appropriate number of points to be circular left
    shifted. However, since then, newer versions of hardware (DQD) and software
    has introduced a new set of values. Depending on the firmware versions
    (DSPFVS) and the decimation rate (DECIM), the following lookup table will
    give the circular shift values needed to correct the DMX data. The values
    of DECIM and DSPFVS can be found in the acqus file in the directory containing
    the data.

     DECIM           DSPFVS 10       DSPFVS 11      DSPFVS 12

       2              44.7500         46.0000        46.311
       3              33.5000         36.5000        36.530
       4              66.6250         48.0000        47.870
       6              59.0833         50.1667        50.229
       8              68.5625         53.2500        53.289
      12              60.3750         69.5000        69.551
      16              69.5313         72.2500        71.600
      24              61.0208         70.1667        70.184
      32              70.0156         72.7500        72.138
      48              61.3438         70.5000        70.528
      64              70.2578         73.0000        72.348
      96              61.5052         70.6667        70.700
     128              70.3789         72.5000        72.524
     192              61.5859         71.3333
     256              70.4395         72.2500
     384              61.6263         71.6667
     512              70.4697         72.1250
     768              61.6465         71.8333
    1024              70.4849         72.0625
    1536              61.6566         71.9167
    2048              70.4924         72.0313


    The number of points obtained from the table are usually not integers.  
    The appropriate procedure is to circular shift (see protocol for details) by the integer
    obtained from truncation of the obtained value and then the residual 1st order phase shift
    that needs to be applied can be obtained by multiplying the decimal
    portion of the calculated number of points by 360.

    For example,

    If DECIM = 32, and DSPFVS = 10,
    then #points 70.0156

    The circular shift performed on the data should be 70 complex points and the linear
    phase correction after Fourier transformation is approximately 0.0156*360 = 5.62 degrees.

    Protocol:

       1. Circular shift (rotate) the appropriate number of points in the data indicated by
       the  DECIM parameter. (see above formulae).

       2. After the circular shift, resize the data to the original size minus
       the number of shifted points. This will leave only the part of the
       data that looks like an FID. Baseline correction (BC) and/or apodization
       (EM etc.) should be applied only on this data, otherwise "In come the frowns."

       Since the first part of the data (the points that are shifted) represents
       negative time, a correct apodization would also multiply the shifted points
       by a negative time apodization. The data size is now returned to
       its original size to reincorporate the shifted points. There may
       still be a discontinuity between the FID portion and the shifted points
       if thelast point of the FID portion is not at zero. This will cause
       sinc wiggles in the peaks.

       3. Applying a zero fill to this data will lead to a discontinuity in the data
       between the rising portion of the shifted points and the zero padding.
       To circumvent this problem, the shifted points are returned (by circular
       shift) to the front of the data, the data is zero filled, and then the
       first points are shifted again to the end of the zero filled data.

       4) The data can now be Fourier transformed and the residual calculated
       1st order phase correction can be applied.
"""
=#
