
using Distributions
using StatsPlots
using Turing
using Logging
using LaTeXStrings

default(labels=false)
     

ways = [0, 3, 8, 9, 0]
ways = ways ./ sum(ways)
println(ways)

b = Binomial(9, 0.5)
pdf(b, 6)

# size of the grid
_bsize = 20

# grid and prior
p_grid = range(0, 1; length=_bsize)
prior = repeat([1.0], _bsize)

# compute likelihood at each value in grid
likelyhood = [pdf(Binomial(9, p), 6) for p in p_grid]

# compute product of likelihood and prior
unstd_posterior = likelyhood .* prior

# standardize the posterior, so it sums to 1
posterior = unstd_posterior / sum(unstd_posterior);


plot(p_grid, posterior; 
    xlabel="probability of water", 
    ylabel="posterior probability",
    title="$size points",
    markershape=:circle
)



@model function water_land(W, L)
    p ~ Uniform(0, 1)
    W ~ Binomial(W + L, p)
end

Logging.disable_logging(Logging.Warn)
chain = sample(water_land(6, 3), NUTS(0.65), 1000)
display(chain)