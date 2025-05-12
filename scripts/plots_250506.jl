# Make some plots from 300 MHz data for May monthly report and posters for 2025
using NMR, GLMakie, CairoMakie
GLMakie.activate!()
import NMR: SFO1, GMAX

datapath = joinpath(homedir(), "nmrdata/data/2025/b300b10/20250430/b300b10")
data = multiwrap(datapath)
#const Dpath = joinpath(homedir(), "nmrdata/data/2025/b800mib/20250423-nug/5")
#dosydata = BrukerSpectrum(Dpath)

# All selective pulses are power-calibrated and have 15 kHz bandwidth.
sel_band = 15e3
# The spectra are in sets of 5: 
# 21-25: GS4 map ReBURP
# 26-30: GS4 ReBURP slice variation, gpz1 0:5:20%
# 31-35: GS4 ReBURP slice variation, no read gradient
# 36-40: as above with cw suppression (residual)
# 41-45: cw on, readgrad on
# 46-49: 15 kHz Gaussian slicer, 5:5:20%

gpz = range(0, 0.2, 5) # rel
gmax= GMAX["BBO-300"].z # G cm⁻¹
sel_cm = get_SliceLength.(gpz, sel_band; Gmax = gmax, message=false)
# prettier version for plotting later
sel_mm = round.(10*sel_cm; sigdigits = 2)
slice_labels = string.(sel_mm) .* " mm, " .* string.(gpz) .* "%"

plt1 = lines(cpu!.(data[21:25]); idxs = (4.7-5., 4.7+5.))

