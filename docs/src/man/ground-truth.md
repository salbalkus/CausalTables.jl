# Computing Ground Truth Conditional Distributions

Once we've defined a DGP and have some table of data with variables matching those of our DGP, we can compute the ground truth conditional distributions of any variable in a CausalTable (given a corresponding DGP) using the `condensity` function. This returns a Distribution object from the package [Distributions.jl](https://juliastats.org/Distributions.jl/stable/).

Let's see an example. First, we'll define a DGP:

```jldoctest truthtest; output = false, filter = r"(?<=.{21}).*"s
using Graphs
using CausalTables
using Random

distributions = @dgp(
        W ~ Binomial(10, 0.3),
        X ~ (@. Normal(:W + 1)),
        Xs = Sum(:X),
        Y ~ (@. LogNormal(log(0.2 * :Xs + 4), 0.1 * :W + 1))
    )

dgp = DataGeneratingProcess(
    n -> erdos_renyi(n, 0.5),
    distributions;
    treatment = :Xs,
    response = :Y,
    controls = [:W]
)

# output
DataGeneratingProcess
```

Now, let's generate some data and compute the ground truth conditional distributions of the variables in the data. Note that if the DGP attempts to summarize a varaible with no neighbors in a graph, the resulting conditional distribution will currently be `Binomial(0, 0.5)`, which denotes a point-mass distribution at 0 in the language of `Distributions.jl`.

```jldoctest truthtest
Random.seed!(1);
data = rand(dgp, 5)
W_distribution = condensity(dgp, data, :W)
X_distribution = condensity(dgp, data, :X)
Xs_distribution = condensity(dgp, data, :Xs)

# output
5-element Vector{Distributions.Normal{Float64}}:
 Distributions.Normal{Float64}(μ=7.0, σ=1.4142135623730951)
 Distributions.Normal{Float64}(μ=8.0, σ=1.4142135623730951)
 Distributions.Normal{Float64}(μ=8.0, σ=1.4142135623730951)
 Distributions.Normal{Float64}(μ=8.0, σ=1.4142135623730951)
 Distributions.Normal{Float64}(μ=9.0, σ=1.4142135623730951)
```

One can also compute the ground truth conditional mean of a variable in a CausalTable using the `conmean` function:

```jldoctest truthtest
Y = conmean(dgp, data, :Y)

# output
5-element Vector{Float64}:
 11.219977823725051
 12.65870031931673
 14.090137140943915
 11.059984096639838
 14.816265013726818
```