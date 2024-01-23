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
using Random

distributions = @dgp(
        W ~ DiscreteUniform(1, 5),
        X ~ (@. Normal(:W, 1)),
        Y ~ (@. Normal(:X + 0.2 * :W, 1))
    )

dgp = DataGeneratingProcess(
    distributions;
    treatment = :X,
    response = :Y,
    controls = [:W]
)

# output
DataGeneratingProcess(CausalTables.var"#26#27"(), Pair{Symbol, Union{NetworkSummary, Function}}[:W => CausalTables.var"#36#38"(), :X => CausalTables.var"#39#41"(), :Y => CausalTables.var"#42#44"()], :X, :Y, [:W])
```

One we've defined our list of distribution functions, we can generate data from the DGP using the `rand` function:

```jldoctest quicktest
Random.seed!(1);
data = rand(dgp, 5)

# output
CausalTable((W = [1, 2, 4, 4, 5], X = [1.267564418628885, 3.749933692528245, 3.1739792080807026, 2.957247582108903, 4.670866154143596], Y = [0.9853125032197694, 5.332176395928403, 4.454670622625683, 3.7546894015545953, 7.105478705857513]), :X, :Y, [:W], SimpleGraph{Int64}(0, Vector{Int64}[]), NamedTuple())
```

For a more detailed guide of how to generate data please refer to [Generating Data](man/generating-data.md).

### Computing Ground Truth Conditional Distributions

Once we've defined a DGP and have some table of data with variables matching those of our DGP, we can compute the ground truth conditional distributions of any variable in a CausalTable (given a corresponding DGP) using the `condensity` function. This returns a Distribution object from the package [Distributions.jl](https://juliastats.org/Distributions.jl/stable/)

```jldoctest quicktest
X_distribution = condensity(dgp, data, :X)

# output
5-element Vector{Normal{Float64}}:
 Normal{Float64}(μ=1.0, σ=1.0)
 Normal{Float64}(μ=2.0, σ=1.0)
 Normal{Float64}(μ=4.0, σ=1.0)
 Normal{Float64}(μ=4.0, σ=1.0)
 Normal{Float64}(μ=5.0, σ=1.0)
```

For convenience, there also exists a `conmean` function that extracts the true conditional mean of a specific variable the CausalTable:

```jldoctest quicktest
Y_mean = conmean(dgp, data, :Y)

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

```jldoctest quicktest
tbl = (W = rand(1:5, 10), X = randn(10), Y = randn(10))
ctbl = CausalTable(tbl; treatment = :X, response = :Y, controls = [:W])

# output
CausalTable((W = [4, 1, 1, 2, 1], X = [0.7979290545612041, -1.2452663162861664, -0.5043635727484704, 0.7380730015204348, -0.5070072144278108], Y = [1.5526026771305754, 0.8274757264628023, 2.121981455933453, -1.5729323029574562, -1.3901156800484877]), :X, :Y, [:W], SimpleGraph{Int64}(0, Vector{Int64}[]), NamedTuple())
```

For a more detailed guide of how to wrap an existing table as a CausalTable please refer to [Turning Your Data Into a CausalTable](man/formatting.md).





