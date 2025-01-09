# Computing Ground Truth Conditional Distributions

One main goal of CausalTables.jl is to allow statisticians to easily develop and test causal inference methods on simulated data. To this end, the package provides a way to compute "ground truth" causal quantities under a particular structural causal model (SCM) and intervention.

Quick methods to approximate the most common ground truth estimands such as the average treatment effect are provided in the section [Approximating ground truth causal estimands](estimands.md). However, in many cases it is also helpful to know the ground truth of more complicated functions of the data, such as conditional distributions, conditional means, or conditional variances. For example, if one is building a machine learning model to predict a conditional mean $\mathbb{E}(Y \mid X = x)$ or a conditional density $p(y \mid X = x)$, it may be useful to know the true conditional mean to evaluate the model's performance. This section will explain how to compute these ground truth quantities using CausalTables.jl.

## Intervening on a `CausalTable`

First, let us set up a simple example. First, we'll define a DGP using the `@dgp` macro and create a StructuralCausalModel object from it. This DGP contains a binary treatment $A$ that depends on a continuous confounder $W$ and a continuous outcome $Y$ that depends on both $A$ and $W$. 

```jldoctest truthtest; output = false, filter = r"(?<=.{21}).*"s
using CausalTables
using Random
using Distributions

# Define the sequence of random variables to be drawn
dgp = @dgp(
    W ~ Beta(2, 4),
    A ~ Bernoulli.(0.5 .* W .+ 0.2),
    Y ~ Normal.(W .+ A, 1)
)

scm = StructuralCausalModel(
    dgp;
    treatment = :A,
    response = :Y
)

# output
StructuralCausalModel
```

In the causal inference setting, we are typically concerned with outcomes under one or more interventions of interest. For example, for binary treatments, we often want to compare the outcome mean had every unit been treated versus had every unit not been treated. CausalTables.jl allows users to easily apply an intervention of interest to a given CausalTable, and subsequently compute counterfactual functions of interest, such as conditional means under intervention.

```jldoctest truthtest; output = false, filter = r"(?<=.{11}).*"s
# Draw a random CausalTable from the StructuralCausalModel
Random.seed!(1);
ct = rand(scm, 500)

# Apply treated and untreated interventions to the CausalTable
ct_treated = intervene(ct, treat_all)
ct_untreated = intervene(ct, treat_none)

# output
CausalTable
```

To evaluate the causal effects of continuous treatments, one can apply the `additive_mtp` and `multiplicative_mtp` functions, which apply an additive or multiplicative shift, respectively, to the natural value of treatment. Their counterfactual differences (compared to the mean response in the source data) provide the causal analogues of a risk difference and risk ratio commonly estimated using linear models in scientific studies. 

## Obtaining Conditional Distributions and Functionals

Once we've defined an SCM (see [Generating data for statistical experiments](generating-data.md)) and have some table of intervened data with variables matching those of our DGP, we can compute the ground truth conditional distributions of any variable (given a corresponding DGP) using the `condensity` function. For any line starting with a `~` in the DGP, `condensity` will be able to return the true conditional distribution, a Distribution object from the package [Distributions.jl](https://juliastats.org/Distributions.jl/stable/), given the data. Some examples are shown below:

```jldoctest truthtest; output = false, filter = r"(?<=.{18}).*"s

# Distribution of the treatment in the observed data
A_distribution = condensity(scm, ct, :A)

# Distribution of the outcome had everyone been treated
Y_under_treatment = condensity(scm, ct_treated, :Y)

# Distribution of the outcome had no one been treated
Y_under_treatment = condensity(scm, ct_untreated, :Y)

# output
500-element Vector
```

One can also compute the ground truth of various functions of these distributions, including the conditional mean (`conmean`), conditional variance (`convar`), or propensity scores (`propensity`; this is the density function evaluated at the observed value of the given variable). To compute other functions of conditional densities not included in CausalTables.jl, please see [Distributions.jl](https://juliastats.org/Distributions.jl/stable/). Below, we show two examples of how one might compute an average treatment effect (ATE) using two different ground-truth functionals of the data.

```jldoctest truthtest; output = false, filter = r"(?<=.{16}).*"s

### Plug-in Estimate ###
μ_treated = conmean(scm, ct_treated, :Y) 
μ_untreated = conmean(scm, ct_untreated, :Y) 
plugin = mean(μ_treated .- μ_untreated) 

### Inverse Propensity Weighted Estimate
p = propensity(scm, ct, :A) 
y = responsematrix(ct)
ipw = mean(y ./ p) .- mean(y) 

# output
1.0037537220251131
```

## Implementing Your Own Interventions

Forthcoming.

## Network Summaries

In causal inference settings featuring *network interference*, the outcome of one unit is not just affected by its own treatment, but also by some summary function of neighboring treatments. CausalTables.jl can also compute conditional densities for certain classes of summary functions. 

If the StructuralCausalModel contains a summary function on a line starting with `$`, `condensity` will return the conditional distribution of the variable being summarized only if it admits a closed-form solution (otherwise an error will be thrown). Note that if the DGP attempts to summarize a variable with no neighbors in a graph, the resulting conditional distribution will currently be `Binomial(0, 0.5)`, which denotes a point-mass distribution at 0. Let's see an example.

```jldoctest truthtest; output = false, filter = r"(?<=.{16}).*"s
using Graphs

dgp = @dgp(
        W ~ Binomial(10, 0.3),
        X ~ (@. Normal(W + 1)),
        A = Graphs.adjacency_matrix(barabasi_albert(length(X), 2)),
        Xs $ Sum(:X, :A),
        Y ~ (@. LogNormal(log(0.2 * Xs + 4), 0.1 * W + 1))
    )
scm = StructuralCausalModel(dgp; treatment = :X, response = :Y)

ct = rand(scm, 5)
W_distribution = condensity(scm, ct, :W)
X_distribution = condensity(scm, ct, :X)
Xs_distribution = condensity(scm, ct, :Xs)

# output
5-element Vector
```

## API

```@autodocs; canonical=false
Modules = [CausalTables]
Order   = [:type, :function]
Pages = ["conditional_density.jl"]
```