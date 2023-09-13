using CSV
using DataFrames
using Turing
using Logging
using StatisticalRethinking 
using StatisticalRethinking: link # import explicitly, because Turing has link method also
using StatisticalRethinkingPlots
using StatsBase
using Random
using LaTeXStrings
using StatsPlots

using Dagitty

using DataFramesMeta, CategoricalArrays



#d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/CC/cubical_entropy.csv"))
d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/Results_fromJulia/GWAS_DataFrame_volume_entropy_Normalized_July2023.csv"))

d = sort!(d, :Volume)
@transform!(d, :DGRP = categorical(:DGRP))
@transform!(d, :sex = categorical(:sex))


d[!,:Volume] = standardize(ZScoreTransform, d.Volume)
d[!,:Entropy0] = standardize(ZScoreTransform, d.Entropy0)
d[!,:Entropy1] = standardize(ZScoreTransform, d.Entropy1)
d[!,:Entropy2] = standardize(ZScoreTransform, d.Entropy2)

d[!,:DGRP_id] = indexin(d.DGRP, unique(d.DGRP));
dgrp_counts = maximum(levels(d.DGRP_id))



@df d scatter(:DGRP_id, :Volume; xlab="DGRP", ylab="Volume (std)",
legend=false, ytickfontsize=8, yticks = :all, size=(800,600))

Random.seed!(100)

@df d scatter3d(:Entropy0, :Entropy1, :Entropy2; xlab="Entropy 0", ylab="Entropy 1", zlab = "Entropy 2", 
camera=(50,20),legend=false, ytickfontsize=8, yticks = :all, size=(800,600), markercolor = "darkred", markersize = 5)
plot!(size=(700,700))

################################################################ Volume model #################################


#d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/Meetings/POLS/Results/new_full_res_June2023_Volume_HRatio.csv"))
d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/Results_fromJulia/GWAS_DataFrame_volume_entropy_Normalized_July2023.csv"))
d = sort!(d, :Volume)
@transform!(d, :DGRP = categorical(:DGRP))
@transform!(d, :sex = categorical(:sex))

d[!,:Volume] = standardize(ZScoreTransform, d.Volume)
d[!,:DGRP_id] = indexin(d.DGRP, unique(d.DGRP));
dgrp_counts = maximum(levels(d.DGRP_id))


@model function model_m5_9(DGRP_id, Volume)
    dgrp_μ = zeros(dgrp_counts)
    a ~ MvNormal(dgrp_μ, 10)
    σ ~ Exponential(1)
    Volume ~ MvNormal(a[DGRP_id], σ)
end

m5_9 = sample(model_m5_9(d.DGRP_id, d.Volume), NUTS(),MCMCThreads(), 1000,4)
plot(m5_9)
m5_9_df = DataFrame(m5_9)
# get rid of square brackets in a params
rename!(n -> replace(n, r"\[|\]" => ""), m5_9_df)

pars = [:a1, :a2, :a3, :a4, :a5, :a6, :a7, :a8, :a9, 
:a10, :a11, :a12, :a13, :a14, :a15, :a16, :a17, :a18, :a19,:a20, :a21, :a22, :a23, 
:a24, :a25, :a26, :a27, :a28, :a29, :a30, :a31, :a32, :a33, :a34, :a34, :a35,
:a36, :a37, :a38, :a39, :a40, :a41, :a42, :a43, :a44, :a45, :a46, :a47, :a48, :a49,
:a50, :a51,:a52, :a53, :a54, :a55, :a56 ] 
p_names = map(v -> "$(v[1]): $(v[2])", zip(pars, unique(d.DGRP)))


df = m5_9_df
coeftab_plot(df; pars=pars, pars_names=p_names, xlab="expected Volume (std)", 
legend=false, ytickfontsize=8, yticks = :all, size=(600,800))


precis(m5_9_df)

post = sample(resetrange(m5_9), 2000)
post = DataFrame(post)
describe(post)




################################################################

d = sort!(d, :Entropy0)
d[!,:DGRP_id] = indexin(d.DGRP, unique(d.DGRP));
dgrp_counts = maximum(levels(d.DGRP_id))



@model function model_m5_10(DGRP_id, Entropy0)
    dgrp_μ = zeros(dgrp_counts)
    a ~ MvNormal(dgrp_μ, 10)
    σ ~ Exponential(1)
    Entropy0 ~ MvNormal(a[DGRP_id], σ)
end

m5_9 = sample(model_m5_10(d.DGRP_id, d.Entropy0), NUTS(),MCMCThreads(), 1000,4)
plot(m5_9)
m5_9_df = DataFrame(m5_9)
# get rid of square brackets in a params
rename!(n -> replace(n, r"\[|\]" => ""), m5_9_df)

pars = [:a1, :a2, :a3, :a4, :a5, :a6, :a7, :a8, :a9, 
:a10, :a11, :a12, :a13, :a14, :a15, :a16, :a17, :a18, :a19,:a20, :a21, :a22, :a23, 
:a24, :a25, :a26, :a27, :a28, :a29, :a30, :a31, :a32, :a33, :a34, :a34, :a35, :a35, :a36, :a37,
:a38, :a39, :a40, :a41, :a42, :a43, :a44, :a45, :a46, :a47, :a48, :a49, :a50, :a51, :a52,
:a53, :a54, :a54, :a55, :a56, :a57, :a58] 
p_names = map(v -> "$(v[1]): $(v[2])", zip(pars, unique(d.DGRP)))


df = m5_9_df
coeftab_plot(df; pars=pars, pars_names=p_names, xlab="expected Entropy0 (std)", 
legend=false, ytickfontsize=8, yticks = :all, size=(600,650))


################################################################
DGRP_seq = range(-2.5, 3.5; length=30)

μ = link(m2_df, [:α, :βE], M_seq)
μ = hcat(μ...)
μ_mean = mean.(eachcol(μ))
μ_PI = PI.(eachcol(μ))
μ_PI = vcat(μ_PI'...)
###################################################################




@model function model_m2(Volume, Entropy0)
    σ ~ Exponential(1)
    α ~ Normal(0, 0.2)
    βE ~ Normal(0, 0.5)
    μ = @. α + βE * Entropy0
    Volume ~ MvNormal(μ, σ)
end 

m2 = sample(model_m2(d.Volume, d.Entropy0), NUTS(),MCMCThreads(), 1000, 4)
plot(m2)
m2_df = DataFrame(m2)

M_seq = range(-2.5, 3.5; length=30)

μ = link(m2_df, [:α, :βE], M_seq)
μ = hcat(μ...)
μ_mean = mean.(eachcol(μ))
μ_PI = PI.(eachcol(μ))
μ_PI = vcat(μ_PI'...)

@df d scatter(:Volume, :Entropy0; xlab="Volume (std)", ylab="Entropy0 (std)", legend=false)
plot!(M_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.2)

coeftab_plot(m2_df)

plot(m2)






########################################################################
# Entropy1 BrainAnalysis


d = sort!(d, :Entropy1)
d[!,:DGRP_id] = indexin(d.DGRP, unique(d.DGRP));
dgrp_counts = maximum(levels(d.DGRP_id))



@model function model_m3(DGRP_id, Entropy1)
    dgrp_μ = zeros(dgrp_counts)
    a ~ MvNormal(dgrp_μ, 10)
    σ ~ Exponential(1)
    Entropy1 ~ MvNormal(a[DGRP_id], σ)
end

m3 = sample(model_m3(d.DGRP_id, d.Entropy1), NUTS(),MCMCThreads(), 1000,4)
plot(m3)
m3_df = DataFrame(m3)
# get rid of square brackets in a params
rename!(n -> replace(n, r"\[|\]" => ""), m3_df)


pars = [:a1, :a2, :a3, :a4, :a5, :a6, :a7, :a8, :a9, 
:a10, :a11, :a12, :a13, :a14, :a15, :a16, :a17, :a18, :a19,:a20, :a21, :a22, :a23, 
:a24, :a25, :a26, :a27, :a28, :a29, :a30, :a31, :a32, :a33, :a34, :a34, :a35, :a35, :a36, :a37,
:a38, :a39, :a40, :a41, :a42, :a43, :a44, :a45, :a46, :a47, :a48, :a49, :a50, :a51, :a52,
:a53, :a54, :a54, :a55, :a56, :a57, :a58] 
p_names = map(v -> "$(v[1]): $(v[2])", zip(pars, unique(d.DGRP)))


df = m3_df
coeftab_plot(df; pars=pars, pars_names=p_names, xlab="expected Entropy1 (std)", 
legend=false, ytickfontsize=8, yticks = :all, size=(600,650))


########################################################################
# Entropy2 BrainAnalysis


d = sort!(d, :Entropy2)
d[!,:DGRP_id] = indexin(d.DGRP, unique(d.DGRP));
dgrp_counts = maximum(levels(d.DGRP_id))



@model function model_m4(DGRP_id, Entropy2)
    dgrp_μ = zeros(dgrp_counts)
    a ~ MvNormal(dgrp_μ, 10)
    σ ~ Exponential(1)
    Entropy2 ~ MvNormal(a[DGRP_id], σ)
end

m4 = sample(model_m4(d.DGRP_id, d.Entropy2), NUTS(),MCMCThreads(), 1000,4)
m4_df = DataFrame(m4)
# get rid of square brackets in a params
rename!(n -> replace(n, r"\[|\]" => ""), m4_df)

pars = [:a1, :a2, :a3, :a4, :a5, :a6, :a7, :a8, :a9, 
:a10, :a11, :a12, :a13, :a14, :a15, :a16, :a17, :a18, :a19,:a20, :a21, :a22, :a23, 
:a24, :a25, :a26, :a27, :a28, :a29, :a30, :a31, :a32, :a33, :a34, :a34, :a35, :a35, :a36, :a37,
:a38, :a39, :a40, :a41, :a42, :a43, :a44, :a45, :a46, :a47, :a48, :a49, :a50, :a51, :a52,
:a53, :a54, :a54, :a55, :a56, :a57, :a58] 
p_names = map(v -> "$(v[1]): $(v[2])", zip(pars, unique(d.DGRP)))


df = m4_df
coeftab_plot(df; pars=pars, pars_names=p_names, xlab="expected Entropy2 (std)", 
legend=false, ytickfontsize=8, yticks = :all, size=(600,650)) 





########################################################################
# StatsPlots


### Organizes the data frame for GWAS`


DGRP_lines_corr = ["DGRP_" * string(i) for i in d.DGRP_line ]
Sex_corr = Vector()

for i in d.Sex
    if 'f' in i 
        push!(Sex_corr,'F')
    else 
        push!(Sex_corr, 'M')
    end
end


dat_gwas = DataFrame(DGRP = DGRP_lines_corr, sex = Sex_corr, Entropy0 = d.Entropy0, Entropy1 = d.Entropy1, Entropy2 = d.Entropy2, Volume = d.Volume)
#CSV.write("/Users/skumar/Documents/PhD/BrainAnalysis/Results_fromJulia/GWAS_DataFrame.csv", dat_gwas)

dat_gwas_volonly = DataFrame(DGRP = DGRP_lines_corr, sex = Sex_corr, Volume = d.Volume)

#CSV.write("/Users/skumar/Documents/PhD/BrainAnalysis/Results_fromJulia/GWAS_DataFrame_Vol_only.csv", dat_gwas_volonly)



########################################################################
# plotting volume and Entropy


@df d scatter3d(:Entropy0, :Entropy2, :Volume; xlab="Entropy 0", ylab="Entropy 2", zlab = "Volume", 
camera=(50,10),legend=false, ytickfontsize=8, yticks = :all, size=(800,600), markercolor = "darkred", markersize = 5)
plot!(size=(700,700))
surface!(d.Entropy0, d.Entropy2, d.Volume)


######################################################################## GWAS only volume

#d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/Results_fromJulia/GWAS_DataFrame_volume_entropy_Normalized_July2023.csv"))
d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/Results_Vol_Entropy/SVM_Entropies_Volume.csv"))



DGRP_lines_corr = ["DGRP_" * string(i) for i in d.DGRP ]
Sex_corr = Vector()

for i in d.sex
    if 'f' in i 
        push!(Sex_corr,'F')
    else 
        push!(Sex_corr, 'M')
    end
end


dat_gwas = DataFrame(DGRP = DGRP_lines_corr, sex = d.sex, 
Volume = d.Volume, Entropy0 = d.Entropy0, Entropy1 = d.Entropy1, Entropy2 = d.Entropy2, Scores = d.scores)
#CSV.write("/Users/skumar/Documents/PhD/BrainAnalysis/Results_fromJulia/GWAS_DataFrame_SVM.csv", dat_gwas)

#dat_gwas_volonly = DataFrame(DGRP = DGRP_lines_corr, sex = Sex_corr, Volume = d.Volume)

#CSV.write("/Users/skumar/Documents/PhD/BrainAnalysis/Results_fromJulia/GWAS_DataFrame_Vol_only.csv", dat_gwas_volonly)



########################################################################

d = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/Results_Vol_Entropy/Results/new_full_res_June2023_Volume_HRatio.csv"))

d2 = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/Results_Vol_Entropy/Results/Entropy_all_July_20_2023_V3_Normalized_Entropy.csv"))

DGRP_lines = Vector()
Sex = Vector()
Entropy0= Vector()
Entropy1 = Vector()
Entropy2 = Vector()
Volume = Vector()


for i in 1:length(d.DGRP_line)
    line = d.DGRP_line[i]
    sex = d.Sex[i]
    print(line)
    for ii in 1:length(d2.DGRP)
        if d2.DGRP[ii] == line && d2.Sex[ii] == sex
            print("Ok")
            push!(DGRP_lines,line)
            push!(Sex,sex)
            push!(Entropy0,d2.entropy0[ii])
            push!(Entropy1,d2.entropy1[ii])
            push!(Entropy2,d2.entropy2[ii])
            push!(Volume,d.Abs_volume[i])
        end
    end
end



dat_Entropy_Vol_full = DataFrame(DGRP = DGRP_lines, Sex = Sex, Volume = Volume, 
Entropy0 = Entropy0, Entropy1 = Entropy1, Entropy2 = Entropy2)

#CSV.write("/Users/skumar/Documents/PhD/BrainAnalysis/Results_fromJulia/GWAS_DataFrame_volume_entropy_Normalized_July2023.csv", dat_Entropy_Vol_full)

