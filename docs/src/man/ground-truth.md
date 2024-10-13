# Computing Ground Truth Conditional Distributions

One main goal of CausalTables.jl is to allow statisticians to easily develop and test causal inference methods on simulated data. To this end, the package provides a way to compute "ground truth" causal quantities under a particular data-generating process (DGP). 

Quick methods to approximate the most common ground truth estimands such as the average treatment effect are provided in the section [Approximating ground truth causal estimands](man/estimands.md). However, in many cases it is also helpful to know the ground-truth of more complicated functions of the data, such as conditional distributions, conditional means, or conditional variances. For example, if one is building a machine learning model to predict a conditional mean $\mathbb{E}(Y \mid X = x)$ or a conditional density $p(y \mid X = x)$, it may be useful to know the true conditional mean to evaluate the model's performance. This section will explain how to compute these ground truth quantities using CausalTables.jl.

Once we've defined a DGP and have some table of data with variables matching those of our DGP, we can compute the ground truth conditional distributions of any variable in a CausalTable (given a corresponding DGP) using the `condensity` function. This returns a Distribution object from the package [Distributions.jl](https://juliastats.org/Distributions.jl/stable/).

Let's see an example. First, we'll define a DGP using the `@dgp` macro and create a StructuralCausalModel object from it. In this DGP, the outcome $Y$ is dependent on the sum of its neighbors' treatments in a random network, meaning that the distributions of some variables in the DGP can be quite complicated. 

```jldoctest truthtest; output = false, filter = r"(?<=.{21}).*"s
using Graphs
using CausalTables
using Random
using Distributions

dgp = @dgp(
        W ~ Binomial(10, 0.3),
        X ~ (@. Normal(W + 1)),
        A = Graphs.adjacency_matrix(barabasi_albert(length(X), 2)),
        Xs $ Sum(:X, :A),
        Y ~ (@. LogNormal(log(0.2 * Xs + 4), 0.1 * W + 1))
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

Now, let's generate some data and compute the ground truth conditional distributions of the variables in the data. For any line starting with a `~` in the DGP, the `condensity` function will be able to return the true conditional distribution of that variable given the data. If the line is a summary function starting with `$`, `condensity` will return the conditional distribution of the variable being summarized only if it admits a closed-form solution (otherwise an error will be thrown). Note that if the DGP attempts to summarize a variable with no neighbors in a graph, the resulting conditional distribution will currently be `Binomial(0, 0.5)`, which denotes a point-mass distribution at 0.

```jldoctest truthtest; output = false, filter = r"(?<=.{16}).*"s
Random.seed!(1);
ctbl = rand(scm, 5)
W_distribution = condensity(scm, ctbl, :W)
X_distribution = condensity(scm, ctbl, :X)
Xs_distribution = condensity(scm, ctbl, :Xs)

# output
5-element Vector
```

One can also compute the ground truth conditional mean or variance of a variable in a CausalTable using the `conmean` and `convar` convenience functions:

```jldoctest truthtest; output = false, filter = r"(?<=.{16}).*"s
σ2 = convar(scm, ctbl, :Y)
μ = conmean(scm, ctbl, :Y)

# output
5-element Vector
```

With these, one can evaluate the performance of flexible statistical models involving functions of conditional densities or moments of a variable in a CausalTable.

### Ground Truth Conditional Distributions API

```@autodocs; canonical=false
Modules = [CausalTables]
Order   = [:type, :function]
Pages = ["conditional_density.jl"]
```