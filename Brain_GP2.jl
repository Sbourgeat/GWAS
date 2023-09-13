using CSV
using DataFrames
using Plots
using GaussianProcesses
using StatisticalRethinking

using DataFramesMeta, CategoricalArrays

using StatsPlots

d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/CC/cubical_entropy.csv"))
d = sort!(d, :Volume)
@transform!(d, :DGRP_line = categorical(:DGRP_line))
@transform!(d, :Sex = categorical(:Sex))


d[!,:Volume] = standardize(ZScoreTransform, d.Volume)
d[!,:Entropy0] = standardize(ZScoreTransform, d.Entropy0)
d[!,:Entropy1] = standardize(ZScoreTransform, d.Entropy1)
d[!,:Entropy2] = standardize(ZScoreTransform, d.Entropy2)

d[!,:DGRP_id] = indexin(d.DGRP_line, levels(d.DGRP_line));

#correlation plot of the dataset
@df d cornerplot([:Entropy0 :Entropy1 :Entropy2 :Volume],grid = false)
plot!(size=(600,600))

# Organize the data
X = Vector()
#push!(X,d[!,:Entropy0])
push!(X,d[!,:Entropy0])
push!(X,d[!,:Entropy2])
Y = d[!,:Volume]


X = mapreduce(permutedims, vcat, X)




mZero = MeanZero()                             # Zero mean function
kern = Matern(5/2,[0.0,0.0],0.0) + SE(0.0,0.0)    # Sum kernel with Matern 5/2 ARD kernel 

gp = GP(X,Y, mZero, kern, -2.0)

optimize!(gp)


dat = vcat([gp.x[1,:]],[gp.x[2,:]],[gp.x[3,:]], [gp.y])
data_gp = DataFrame(dat, [:Entropy0, :Entropy1, :Entropy2, :Volume])


surface(gp.x[1,:],gp.x[2,:],gp.y, 
zcolor = gp.y,xlim=[-2,1.5],zlim=[-5,5],ylim=[-2,2], camera=(30,30))
plot!(xlab="Entropy0", ylab="Entropy2",zlab=("Volume"), colorbar_title="Volume")
plot!(size=(700,700))

plot(contour(gp), heatmap(gp))
plot!(size=(700,700))


wireframe(gp.x[1,:],gp.x[2,:],gp.x[3,:])

 


#Sanity checks

using Distributions

set_priors!(kern, [Normal(), Normal()]) # Uniform(0,1) distribution assumed by default if priors not specified
chain = mcmc(gp)
plot(chain', label=["Noise", "SE log length", "SE log scale"]; fmt=:png)
plot!(size=(800,600))



# predict 

μ, σ² = predict_y(gp, [0,1,2])