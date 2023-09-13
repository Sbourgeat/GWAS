1+1 

exp(22)

x = 2

x < 0 && print("x is negative")

x < 0 || print("x is negative")

function times_2(x)
    return x * 2 
end

times_2(10)


# create a function names compose 
function compose(x, y=10; a, b=10)
    return x, y, a, b
end

compose(1, 2; a=3, b= 4)


using CSV
using DataFrames

df = DataFrame(CSV.File("/Users/skumar/Documents/PhD/BrainAnalysis/CC/cubical_entropy.csv"))

using PlotlyJS


PlotlyJS.plot(df,x=:Entropy0, y=:Entropy1, color=:Sex, mode="markers")





function fun6(x) 
    for i in [1, 2, 3]
        if i == 1
            x = 1
        else 
            x += 1
        end
        println(x)
    end
end




aq = [10.0   8.04  10.0  9.14  10.0   7.46   8.0   6.58
              8.0   6.95   8.0  8.14   8.0   6.77   8.0   5.76
             13.0   7.58  13.0  8.74  13.0  12.74   8.0   7.71
              9.0   8.81   9.0  8.77   9.0   7.11   8.0   8.84
             11.0   8.33  11.0  9.26  11.0   7.81   8.0   8.47
             14.0   9.96  14.0  8.1   14.0   8.84   8.0   7.04
              6.0   7.24   6.0  6.13   6.0   6.08   8.0   5.25
              4.0   4.26   4.0  3.1    4.0   5.39  19.0  12.50
             12.0  10.84  12.0  9.13  12.0   8.15   8.0   5.56
              7.0   4.82   7.0  7.26   7.0   6.42   8.0   7.91
              5.0   5.68   5.0  4.74   5.0   5.73   8.0   6.89]



size(aq)




using Statistics

mean(aq; dims=2)




[mean(col) for col in eachcol(aq)]



y = aq[:,2]

X = [ones(11) aq[:,1]]

plot(X[:,2],y,seriestype=:scatter)



function R2(x, y)
    X = [ones(11) x]
               model = X \ y
               prediction = X * model
               error = y - prediction
               SS_res = sum(v -> v ^ 2, error)
               mean_y = mean(y)
               SS_tot = sum(v -> (v - mean_y) ^ 2, y)
               return 1 - SS_res / SS_tot
           end

plot(scatter(aq[:, 1], aq[:, 2]; legend=false),
scatter(aq[:, 3], aq[:, 4]; legend=false),
scatter(aq[:, 5], aq[:, 6]; legend=false),
scatter(aq[:, 7], aq[:, 8]; legend=false))
