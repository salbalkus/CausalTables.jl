# Computing Ground Truth Conditional Distributions

Once we've defined a DGP and have some table of data with variables matching those of our DGP, we can compute the ground truth conditional distributions of any variable in a CausalTable (given a corresponding DGP) using the `condensity` function. This returns a Distribution object from the package [Distributions.jl](https://juliastats.org/Distributions.jl/stable/).

Let's see an example. First, we'll define a DGP:

```jldoctest truthtest; output = false, filter = r"(?<=.{21}).*"s
using Graphs
using CausalTables
using Random

dgp = @dgp(
        W ~ Binomial(10, 0.3),
        X ~ (@. Normal(:W + 1)),
        A = adjacency_matrix(barabasi_albert(length(:X), 2)),
        Xs $ Sum(:X, :A),
        Y ~ (@. LogNormal(log(0.2 * :Xs + 4), 0.1 * :W + 1))
    )

scm = StructuralCausalModel(
    dgp;
    treatment = :X,
    response = :Y,
    confounders = [:W]
)

# output
StructuralCausalModel
```

Now, let's generate some data and compute the ground truth conditional distributions of the variables in the data. Note that if the DGP attempts to summarize a variable with no neighbors in a graph, the resulting conditional distribution will currently be `Binomial(0, 0.5)`, which denotes a point-mass distribution at 0.

```jldoctest truthtest; output = false, filter = r"(?<=.{16}).*"s
Random.seed!(1);
data = rand(scm, 5)
W_distribution = condensity(scm, data, :W)
X_distribution = condensity(scm, data, :X)
Xs_distribution = condensity(scm, data, :Xs)

# output
5-element Vector
```

One can also compute the ground truth conditional mean of a variable in a CausalTable using the `conmean` function:

```jldoctest truthtest; output = false, filter = r"(?<=.{16}).*"s
Y = conmean(scm, data, :Y)

# output
5-element Vector
```