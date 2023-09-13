using CSV
using DataFrames
using Plots
using GaussianProcesses
using StatisticalRethinking

using DataFramesMeta, CategoricalArrays

using StatsPlots

#d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/CC/cubical_entropy.csv"))
#d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/Results_Vol_Entropy/SVM_Entropies_plus_Volume.csv"))
d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/Results_Vol_Entropy/GWAS_Normalized_EntropyVol.tsv"))



d = sort!(d, :Volume)
@transform!(d, :DGRP = categorical(:DGRP))
@transform!(d, :sex = categorical(:sex))


d[!,:Volume] = standardize(ZScoreTransform, d.Volume)
d[!,:Entropy0] = standardize(ZScoreTransform, d.Entropy0)
d[!,:Entropy1] = standardize(ZScoreTransform, d.Entropy1)
d[!,:Entropy2] = standardize(ZScoreTransform, d.Entropy2)


d[!,:Entropy0] = [2^i for i in d.Entropy0]
d[!,:Entropy1] = [2^i for i in d.Entropy1]
d[!,:Entropy2] = [2^i for i in d.Entropy2]



d[!,:DGRP_id] = indexin(d.DGRP, (d.DGRP));



@df d cornerplot([:Entropy0 :Entropy1 :Entropy2 :Volume],grid = false)
plot!(size=(600,600))

# Organize the data
X = Vector()
#push!(X,d[!,:scores])
push!(X,d[!,:Entropy1])
push!(X,d[!,:Entropy2])
#push!(X,d[!,:Entropy2])
Y = d[!,:Volume]


X = mapreduce(permutedims, vcat, X)




mZero = MeanZero()                 # Zero mean function
kern = Matern(5/2,[0.0, 0.0],0.0)   # Sum kernel with Matern 5/2 ARD kernel 

gp = GP(X,Y, mZero, kern, 10.0)

optimize!(gp)

plot(gp)
#wireframe(gp)
#surface(gp)

gr() # use GR backend to allow wireframe plot
p1 = contour(gp)
p2 = surface(gp)
p3 = heatmap(gp)
p4 = wireframe(gp)
plot(p1, p2, p3, p4; fmt=:png)
plot!(size=(800,600))

#dat = vcat([gp.x[1,:]],[gp.x[2,:]],[gp.x[3,:]], [gp.y])
#data_gp = DataFrame(dat, [:Entropy0, :Entropy1, :Entropy2, :Volume])


#Plots.scatter3d(gp.x[1,:],gp.x[2,:],gp.x[3,:], 
#zcolor = gp.y,xlim=[-3,3],zlim=[-10,16],ylim=[-2,2],camera=(50,30))
#plot!(xlab="Entropy0", ylab="Entropy1",zlab=("Entropy2"), colorbar_title="Volume")
#plot!(size=(700,700))



#wireframe(gp.x[1,:],gp.x[2,:],gp.x[3,:])

 


#Sanity checks

using Distributions

set_priors!(kern,[Normal(), Normal()])
chain = mcmc(gp)
plot(chain', label=["Noise", "SE log length", "SE log scale"]; fmt=:png)
plot!(size=(800,600))



# Prediction
X = Vector()
push!(X,d[!,:Entropy0])
push!(X,d[!,:Entropy1])
push!(X,d[!,:Entropy2])


X = mapreduce(permutedims, vcat, X)


res_gp = predict_f(gp,X)

new_d_for_gwas = Vector()
push!(new_d_for_gwas, d.DGRP_line)
push!(new_d_for_gwas, res_gp[1])
X_gwas = mapreduce(permutedims, vcat, new_d_for_gwas)


new_df_for_gwas = DataFrame(DGRP_lines = d.DGRP_line,Sex = d.Sex, Volume = res_gp[1])
#CSV.write("/Users/skumar/Documents/PhD/BrainAnalysis/Results_fromJulia/GaussianRes.csv", new_df_for_gwas)



using Gaston

surf(data_gp.Entropy0, data_gp.Entropy1, data_gp.Entropy2, lc =:turquoise)







d[!, :entropy]