using CSV
using DataFrames
using PlotlyJS
using DataFramesMeta, CategoricalArrays
using StatisticalRethinking 

d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/CC/cubical_entropy.csv"))
d = sort!(d, [:Volume, :Sex])
@transform!(d, :DGRP_line = categorical(:DGRP_line))
@transform!(d, :Sex = categorical(:Sex))



d[!,:Volume] = standardize(ZScoreTransform, d.Volume)
d[!,:Entropy0] = standardize(ZScoreTransform, d.Entropy0)
d[!,:Entropy1] = standardize(ZScoreTransform, d.Entropy1)
d[!,:Entropy2] = standardize(ZScoreTransform, d.Entropy2)



plot(scatter(d, x=:DGRP_line, y=:Volume, group=:Sex,mode = "markers"))


plot(
    d,
    kind ="scatter",
    mode="markers",
    x=:DGRP_line,
    y=:Volume,
    group=:Sex
)


trace = scatter(d, x=:DGRP_line, y=:Volume, group=:Sex
, mode="markers",)
layout = Layout(xaxis_type="category"
)

plot(trace, layout)


trace = scatter3d(d, x=:Entropy0, y=:Entropy1, z=:Volume,group=:Sex
, mode="markers",)


plot(trace)



# INPUT DATA

cols = [:Volume, :Entropy0, :Entropy1, :Entropy2]  # define subset
M = cor(Matrix(d[!,cols]))       # correlation matrix

# PLOT
(n,m) = size(M)
Plots.heatmap(M, fc=cgrad([:white,:dodgerblue4]), xticks=(1:m,cols), xrot=90, yticks=(1:m,cols), yflip=true)
annotate!([(j, i, text(round(M[i,j],digits=3), 8,"Computer Modern",:black)) for i in 1:n for j in 1:m])