# CausalTables.jl

*A package for storing and simulating data for causal inference in Julia.*

## Installation
CausalTables.jl can be installed using the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```
Pkg> add CausalTables
```

## Quick Start

CausalTables.jl has three main functionalities:

1. Generating simulation data using a `StructuralCausalModel`.
2. Computing "ground truth" conditional distributions, moments, counterfactuals, and counterfactual functionals from a `StructuralCausalModel` and a `CausalTable`. These include, for instance, counterfactual means and average treatment effects.
3. Wrapping an existing Table as a `CausalTable` object for use by external packages, which provides several utility functions for extracting causal-relevant variables from a dataset. 

The examples below illustrate each of these three functionalities.

### Simulating Data from a DataGeneratingProcess

To set up a statistical simulation using CausalTables.jl, we first define a `StructuralCausalModel` (SCM). This consists of two parts: a `DataGeneratingProcess` (DGP) that controls how the data is generated, and a list of variables to define the basic structure of the underlying causal diagram.

A DataGeneratingProcess can be constructed using the `@dgp` macro, which takes a sequence of conditional distributions of the form `[name] ~ Distribution(args...)` or auxiliary variables `[name] = some code...` and returns a `DataGeneratingProcess` object. Then, one can construct an StructuralCausalModel by passing the DGP to its construct, along with labels of the treatment and response variables. Note that `using Distributions` is almost always required before defining a DGP, since the package [Distributions.jl](https://juliastats.org/Distributions.jl/stable/) is used to define the conditional distribution of random components at each step.

```jldoctest quicktest; output = false, filter = r"(?<=.{21}).*"s
using CausalTables
using Random
using Distributions

dgp = @dgp(
        W ~ DiscreteUniform(1, 5),
        X ~ (@. Normal(W, 1)),
        Y ~ (@. Normal(X + 0.2 * W, 1))
    )

scm = StructuralCausalModel(
    dgp;
    treatment = :X,
    response = :Y
)

# output
StructuralCausalModel
```

One we've defined our list of distribution functions, we can generate data from the DGP using the `rand` function:

```jldoctest quicktest; output = false, filter = r"(?<=.{11}).*"s
Random.seed!(1);
data = rand(scm, 5)

# output
CausalTable
```

We can also apply various causal interventions to the data using the `intervene` function. The example below computes a new version of the CausalTable with each unit's treatment shifted by 1 -- this is analogous to the effect estimated by a classical linear regression analysis. 

```jldoctest quicktest; output = false, filter = r"(?<=.{11}).*"s
data_intervened = intervene(data, additive_mtp(1))

# output
CausalTable
```

For a more detailed guide of how to generate data please refer to [Generating Data](man/generating-data.md).

### Computing Ground Truth Functionals

Once we've defined a DGP, we can approximate ground truth statistical functionals along with their efficiency bounds (variance of the counterfactual outcome) for a specified SCM using built-in functions. In general, these include

- Counterfactual Means (`cfmean`)
- Counterfactual Differences (`cfdiff`)
- Average Treatment Effect (`ate`), including among the Treated (`att`) and Untreated (`atu`)
- Average Policy Effect (`ape`), also known as the causal effect of a Modified Treatment Policy. 

For the complete API of available ground truth causal estimands, see [Estimands](man/estimands.md)

```jldoctest quicktest
cfmean(scm, additive_mtp(1))

# output
(μ = 4.599337273915866, eff_bound = 4.881412474779794)
```

For problems that involving functionals not available through CausalTables.jl or that require more fine-grained knowledge of the true conditional distributions for a given dataset, this package also implements the `condensity` function. This function computes the true conditional distributions of any variable in a CausalTable (given a corresponding DGP). The function returns a vector of Distribution objects from the package [Distributions.jl](https://juliastats.org/Distributions.jl/stable/)

```jldoctest quicktest
X_distribution = condensity(scm, data, :X)

# output
5-element Vector{Normal{Float64}}:
 Distributions.Normal{Float64}(μ=1.0, σ=1.0)
 Distributions.Normal{Float64}(μ=2.0, σ=1.0)
 Distributions.Normal{Float64}(μ=4.0, σ=1.0)
 Distributions.Normal{Float64}(μ=4.0, σ=1.0)
 Distributions.Normal{Float64}(μ=5.0, σ=1.0)
```

For convenience, there also exists functins like `conmean`, `convar`, and `propensity` that extract the true conditional mean, variance, and (generalized) propensity score of a specific variable the CausalTable. One can apply this to an "intervened" version of data to obtain functionals of the outcome under intervention:

```jldoctest quicktest
Y_mean = conmean(scm, data_intervened, :Y)

# output
5-element Vector{Float64}:
 2.467564418628885
 5.149933692528245
 4.973979208080702
 4.757247582108903
 6.670866154143596
```

For a more detailed guide of how to compute ground truth conditional distributions please refer to [Computing Ground Truth Conditional Distributions](man/ground-truth.md).

### Wrapping an existing Table as a CausalTable

If you have a table of data that you would like to use with CausalTables.jl without defining a corresponding DataGeneratingProcess (i.e. to use with another package) you can wrap it as a `CausalTable` using the corresponding constructor.

```jldoctest quicktest; output = false, filter = r"(?<=.{11}).*"s
tbl = (W = rand(1:5, 10), X = randn(10), Y = randn(10))
ctbl = CausalTable(tbl; treatment = :X, response = :Y, 
                        causes = (X = [:W], Y = [:W, :X]))

# output
CausalTable
```

Doing this is often convenient, as it allows you to use the utility functions provided by CausalTables.jl to extract causal-relevant variables from the dataset. For instance, you can extract the treatment, response, confounders, mediators, or instruments from the dataset using the corresponding functions. For example, the following subsets the data to include only confounders:

```jldoctest quicktest; output = false, filter = r"(?<=.{11}).*"s
confounders(ctbl)

# output
CausalTable
```

For a more detailed guide of how to wrap an existing table as a CausalTable please refer to [Turning Your Data Into a CausalTable](man/formatting.md).

# Contributing

Have questions? Spot a bug or issue in the documentation? Want to request a new feature or add one yourself? Please don't hesitate to open an issue or pull request on the [CausalTables.jl GitHub repository](https://github.com/salbalkus/CausalTables.jl). We welcome all contributions and feedback!
