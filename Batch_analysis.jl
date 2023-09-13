using CSV
using DataFramesMeta, CategoricalArrays


d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/Meetings/POLS/Results/new_full_res_June2023_Volume_HRatio.csv"))
d_batch = DataFrame(CSV.File("/Users/skumar/Documents/PhD/Meetings/POLS/Results/Batch_res.csv"))


@transform!(d, :DGRP_line = categorical(:DGRP_line))
@transform!(d_batch, :DGRP = categorical(:DGRP))


Batch = Vector()
Lines = Vector()
Vol = Vector()
Sex = Vector()
for i in 1:length(d.DGRP_line)
    for ii in 1:length(d_batch.DGRP)
        if d.DGRP_line[i] == d_batch.DGRP[ii]
            push!(Batch, d_batch.Batch[ii])
            push!(Lines, d.DGRP_line[i])
            push!(Vol, d.Abs_volume[i])
            push!(Sex, d.Sex[i])
            
        end
    end
end


data_batch = DataFrame(Lines = Lines, Batch = Batch, Vol = Vol, Sex = Sex)


#CSV.write("/Users/skumar/Documents/PhD/BrainAnalysis/Results_fromJulia/Batch_DataFrame.csv", data_batch)
using Plots
using StatsPlots

@df data_batch Plots.scatter(:Lines,:Vol, group=:Batch, legend=false, size=(800,600))


@df data_batch Plots.boxplot(:Batch,:Vol)


@df data_batch density(:Vol, group =:Batch, legend = :topleft)


data_batch = sort!(data_batch, :Vol)
data_batch[!,:DGRP_id] = indexin(data_batch.Lines, unique(data_batch.Lines));
@df data_batch Plots.scatter(:DGRP_id,:Vol, group=:Batch, legend=false, size=(800,600))


