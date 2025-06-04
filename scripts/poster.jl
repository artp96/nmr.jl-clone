using NMR_acp
using CairoMakie

figpath = "/home/w25612ap/PhD/Poster 2025/portrait/figures/"

path800 = "/home/w25612ap/nmrdata/data/2025/b800mib/20250602/"
path300 = "/home/w25612ap/nmrdata/data/2025/b300b10/20250430/b300b10/"

data800 = multiwrap(path800)
data300 = multiwrap(path300)

ppmH2O = (4.66, 4.74)
noe_v_gs = [13, 17, 18, 19]
labels = ["NOE-pr" "0% Sel. Grd." "1% Sel. Grd." "2% Sel. Grd."]
H2O_800 = lines([data800[s] |> cpu! for s in noe_v_gs]; idxs = ppmH2O,
                     title = "H₂O Residual",
                     labels = labels)

ppmTSP = (0.1, -0.1)
TSP_800 = lines([data800[s] |> cpu! for s in noe_v_gs]; idxs = ppmTSP,
                     title = "0.5 mM TSP Sig.",
                     labels = labels)
