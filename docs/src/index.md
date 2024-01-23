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

1. Generating simulation data using a `DataGeneratingProcess`
2. Computing "ground truth" conditional distributions, means, and other functionals from a `DataGeneratingProcess` and a `CausalTable`
3. Wrapping an existing Table to make it a `CausalTable` for use by external packages.

The examples below illustrate each of these three functionalities.

### Simulating Data from a DataGeneratingProcess

To set up a statistical simulation using CausalTables.jl, we first define a `DataGeneratingProcess` (DGP). The easiest way to do this is using the `@dgp` macro, which takes a sequence of conditional distributions of the form `[variable name] ~ Distribution(args...)` and returns a `DataGeneratingProcess` object like so:

```jldoctest quicktest
using CausalTables

distributions = @dgp(
        W ~ DiscreteUniform(1, 5),
        X ~ (@. Normal(:W, 1)),
        Y ~ (@. Normal(:A + 0.2 * :W, 1))
    )

dgp = DataGeneratingProcess(
    distributions;
    treatment = :X,
    response = :Y,
    controls = [:W]
)
```

One we've defined our list of distribution functions, we can generate data from the DGP using the `rand` function:

```jldoctest quicktest
data = rand(dgp, 10)
```

For a more detailed guide of how to generate data please refer to [Generating Data](man/generating-data.md).

### Computing Ground Truth Conditional Distributions

Once we've defined a DGP and have some table of data with variables matching those of our DGP, we can compute the ground truth conditional distributions of any variable in a CausalTable (given a corresponding DGP) using the `condensity` function. This returns a Distribution object from the package [Distributions.jl](https://juliastats.org/Distributions.jl/stable/)

```jldoctest quicktest
Y_distribution = condensity(dgp, data, :Y)
```

For convenience, there also exists a `conmean` function that extracts the true conditional mean of a specific variable the CausalTable:

```jldoctest quicktest
Z_mean = conmean(dgp, data, :Z)
```

For a more detailed guide of how to compute ground truth conditional distributions please refer to [Computing Ground Truth Conditional Distributions](man/ground-truth.md).

### Wrapping an existing Table as a CausalTable

If you have a table of data that you would like to use with CausalTables.jl without defining a corresponding DataGeneratingProcess (i.e. to use with another package) you can wrap it as a `CausalTable` using the corresponding constructor:

```jldoctest quicktest
tbl = (W = rand(1:5, 10), X = randn(10), Y = randn(10))
ctbl = CausalTable(tbl; treatment = :X, response = :Y, controls = [:W])
```

For a more detailed guide of how to wrap an existing table as a CausalTable please refer to [Turning Your Data Into a CausalTable](man/formatting.md).





