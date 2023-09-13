using Distributions
using StatsPlots
using StatsBase
using LaTeXStrings
using CSV
using DataFrames
using StatisticalRethinking
using StatisticalRethinking: link
using LinearAlgebra
using Logging
using Random
using Turing

################################################################
# If the packages are not installed
import Pkg;
Pkg.add(["Distributions", "StatsPlots", 
"LaTeXStrings", "StatsBase","CSV","DataFrames", 
"StatisticalRethinking", "LinearAlgebra", 
"Logging", "Random", "Turing", "StatisticalRethinkingPlots"])
################################################################


n = rand(Uniform(-1, 1), 1000, 16);
pos = sum.(eachrow(n));
density(pos)


prod(1 .+ rand(Uniform(0, .1), 12))





u = Uniform(0, .1)
growth = prod.(eachrow(1 .+ rand(u, 10000, 12)));

density(growth; label="density")
# overlay normal distribution
μ = mean(growth)
σ = std(growth)
plot!(Normal(μ, σ); label="Normal")
     
big = prod.(eachrow(1 .+ rand(Uniform(0, 0.5), 10000, 12)));
small = prod.(eachrow(1 .+ rand(Uniform(0, 0.01), 10000, 12)));
density([big, small]; layout=(2, 1))

