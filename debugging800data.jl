using NMR, GLMakie, CairoMakie
GLMakie.activate!()
import NMR: SFO1, GMAX

const dpath = joinpath(homedir(), "nmrdata/data/2025/b800mib/20250417")

# load the spectra from that day
data = multiimport(dpath)
expnos = [d.expno for d in data]
S = (; (Symbol(expnos) .=> data)... )

# a function to find the differences between all the stored metadata of two spectra
function getdiffs(s1::D, s2::D) where D <: AbstractDict
    acqudiffs = Vector{Any}[]
    for key in keys(merge(s1, s2))
           if all(s1[key] .== s2[key])
           #println("$key : ok!")
           else
            println("$key : $(acqu[key][1:30]) \n ≠ \n $(acqu[key][1:30])")
            push!(acqudiffs, [key])
       end
    end

    return acqudiffs
end
getdiffs(s1::S, s2::S) where S <: AbstractSpectrum = getdiffs(s1.acqu, s2.acqu)

s30, s31 = data[30:31]


diffs = getdiffs(s30,s31)

fidscale = maximum(abs, s30.fid) / maximum(abs, s31.fid)
prcscale = maximum(abs, s30[:]) / maximum(abs, s31[:])

idx1, idx2 = ppmtoindex(s30, (24., 25.))
baselinescale = norm(s30[idx1:idx2] ./ s31[idx1:idx2])


