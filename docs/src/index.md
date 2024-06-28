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

1. Generating simulation data using a `StructuralCausalModel`
2. Computing "ground truth" conditional distributions, means, and other functionals from a `DataGeneratingProcess` and a `CausalTable`
3. Wrapping an existing Table to make it a `CausalTable` for use by external packages.

The examples below illustrate each of these three functionalities.

### Simulating Data from a DataGeneratingProcess

To set up a statistical simulation using CausalTables.jl, we first define a `StructuralCausalModel` (SCM). This consists of two parts: a `DataGeneratingProcess` (DGP) that controls how the data is generated, and a list of variables to define the basic structure of the underlying causal diagram.

A DataGeneratingProcess can be constructed using the `@dgp` macro, which takes a sequence of conditional distributions of the form `[variable name] ~ Distribution(args...)` and returns a `DataGeneratingProcess` object. Then, one can construct an StructuralCausalModel by passing the DGP to its construct, along with labels of the treatment and response variables.

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
    response = :Y,
    confounders = [:W]
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

For a more detailed guide of how to generate data please refer to [Generating Data](man/generating-data.md).

### Computing Ground Truth Conditional Distributions

Once we've defined a DGP and have some table of data with variables matching those of our DGP, we can compute the ground truth conditional distributions of any variable in a CausalTable (given a corresponding DGP) using the `condensity` function. This returns a Distribution object from the package [Distributions.jl](https://juliastats.org/Distributions.jl/stable/)

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

For convenience, there also exists a `conmean` function that extracts the true conditional mean of a specific variable the CausalTable:

```jldoctest quicktest
Y_mean = conmean(scm, data, :Y)

# output
5-element Vector{Float64}:
 1.467564418628885
 4.149933692528245
 3.973979208080703
 3.757247582108903
 5.670866154143596
```

For a more detailed guide of how to compute ground truth conditional distributions please refer to [Computing Ground Truth Conditional Distributions](man/ground-truth.md).

### Wrapping an existing Table as a CausalTable

If you have a table of data that you would like to use with CausalTables.jl without defining a corresponding DataGeneratingProcess (i.e. to use with another package) you can wrap it as a `CausalTable` using the corresponding constructor:

```jldoctest quicktest; output = false, filter = r"(?<=.{11}).*"s
tbl = (W = rand(1:5, 10), X = randn(10), Y = randn(10))
ctbl = CausalTable(tbl; treatment = :X, response = :Y, confounders = [:W])

# output
CausalTable
```

For a more detailed guide of how to wrap an existing table as a CausalTable please refer to [Turning Your Data Into a CausalTable](man/formatting.md).





