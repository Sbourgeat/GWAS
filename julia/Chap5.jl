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


# Spurious association

## Standardize the data with ZScoreTransform
d = DataFrame(CSV.File("WafflesDivorce.csv"))
d[!,:D] = standardize(ZScoreTransform, d.Divorce)
d[!,:M] = standardize(ZScoreTransform, d.Marriage)
d[!,:A] = standardize(ZScoreTransform, d.MedianAgeMarriage)

std(d.MedianAgeMarriage)

Random.seed!(100)

# @ signifies a macro operator
@model function model_m5_1(A,D)
    σ ~ Exponential(1) 
    α ~ Normal(0, 0.2)
    βA ~ Normal(0, 0.5)
    μ = @. α + βA * A
    D ~ MvNormal(μ, σ)
end

m5_1 = sample(model_m5_1(d.A, d.D), NUTS(), 1000)
m5_1_df = DataFrame(m5_1)
prior = sample(model_m5_1([0],[0]), Prior(), 1000)
prior_df = DataFrame(prior)


# calculate μ for every prior sample on age = -2 and age = 2
bounds = [-2, 2]
μ = link(prior_df, [:α, :βA], bounds)
μ = hcat(μ...);

p = plot(xlab="MedianAgeMarriage (std)", ylab="Divorce rate (std)", legend=false)
for μₚ ∈ first(eachrow(μ), 50)
    plot!(bounds, μₚ; c=:black, alpha=0.3)
end
display(p)       


# Now, plot the linear regression from the posterior distribution

A_seq = range(-3, 3.2; length=30)
μ = link(m5_1_df, [:α, :βA], A_seq)
μ = hcat(μ...)
μ_mean = mean.(eachcol(μ))
μ_PI = PI.(eachcol(μ))
μ_PI = vcat(μ_PI'...)


@df d scatter(:A, :D; xlab="Median age marriage", 
ylab="Divorce rate", legend=false)
plot!(A_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.2)


# Model Mariage rate vs Divorce rate this time

Random.seed!(100)

@model function model_m5_2(M, D)
    σ ~ Exponential(1)
    α ~ Normal(0, 0.2)
    βM ~ Normal(0, 0.5)
    μ = @. α + βM * M
    D ~ MvNormal(μ, σ)
end 

m5_2 = sample(model_m5_2(d.M, d.D), NUTS(), 1000)
m5_2_df = DataFrame(m5_2)

M_seq = range(-1.74, 2.8; length=30)

μ = link(m5_2_df, [:α, :βM], M_seq)
μ = hcat(μ...)
μ_mean = mean.(eachcol(μ))
μ_PI = PI.(eachcol(μ))
μ_PI = vcat(μ_PI'...)

@df d scatter(:M, :D; xlab="Marriage rate (std)", ylab="Divorce rate (std)", legend=false)
plot!(M_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.2)


g = Dagitty.DAG(:A => :M, :A => :D, :M => :D)
drawdag(g, [0, 1, 2], [0, 1, 0])

g = Dagitty.DAG(:A => :M, :A => :D)
implied_conditional_independencies(g)


g = Dagitty.DAG(:A => :M, :A => :D, :M => :D)
implied_conditional_independencies(g)


@model function model_m5_3(A, M, D)
    σ ~ Exponential(1)
    α ~ Normal(0, 0.2)
    βA ~ Normal(0, 0.5)
    βM ~ Normal(0, 0.5)
    μ = @. α + βA * A + βM * M
    D ~ MvNormal(μ, σ)
end

m5_3 = sample(model_m5_3(d.A, d.M, d.D), NUTS(), 1000)
m5_3_df = DataFrame(m5_3)
precis(m5_3_df)
coeftab_plot(m5_3_df)    
coeftab_plot(m5_1_df, m5_2_df, m5_3_df; pars=(:βA, :βM), names=["m5.1", "m5.2", "m5.3"]) 

# Interpretation, once we know about the age, 
# we don't learn anything more looking at 
# the mariage rate concerning Divorce



Random.seed!(100)
N = 50 
age = rand(Normal(), N)
mar = rand.(Normal.(-age))
div = rand.(Normal.(age));

s1 = DataFrame(sample(model_m5_1(age, div), NUTS(), 1000))
s2 = DataFrame(sample(model_m5_2(mar, div), NUTS(), 1000))
s3 = DataFrame(sample(model_m5_3(age,mar, div), NUTS(), 1000));
coeftab_plot(s1, s2, s3; pars=(:βA, :βM) , names=["s1", "s2", "s3"])





Random.seed!(100)
N = 50 
age = rand(Normal(), N)
mar = rand.(Normal.(-age))
div = rand.(Normal.(age .+ mar));

s1 = DataFrame(sample(model_m5_1(age, div), NUTS(), 1000))
s2 = DataFrame(sample(model_m5_2(mar, div), NUTS(), 1000))
s3 = DataFrame(sample(model_m5_3(age,mar, div), NUTS(), 1000));
coeftab_plot(s1, s2, s3; pars=(:βA, :βM) , names=["s1", "s2", "s3"])



# Interaction model

@model function model_m5_4(A, M)
    σ ~ Exponential(1)
    α ~ Normal(0, 0.2)
    βAM ~ Normal(0, 0.5)
    μ = @. α + βAM * A
    M ~ MvNormal(μ, σ)
end

m5_4 = sample(model_m5_4(d.A, d.M), NUTS(), 1000)
m5_4_df = DataFrame(m5_4)

μ = link(m5_4_df, [:α, :βAM], d.A)
μ = hcat(μ...)
μ_mean = mean.(eachcol(μ))
μ_resid = μ_mean .- d.M; 

# Side-note: how to plot the residuals
# getting yerr - list of 2-tuples with distance to the regression line

yerr = collect(zip(-clamp.(μ_resid, -Inf, -0.0), clamp.(μ_resid, 0, Inf)));

plot(d.A, μ_mean, xlab="Age at marriage (std)", ylab="Marriage rate (std)")
scatter!(d.A, d.M)
scatter!(d.A, d.M, yerr=yerr, markersize=0)

