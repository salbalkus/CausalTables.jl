# CausalTables.jl

[![Build Status](https://github.com/salbalkus/CausalTables.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/salbalkus/CausalTables.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/salbalkus/CausalTables.jl/badge.svg?branch=main)](https://coveralls.io/github/salbalkus/CausalTables.jl?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Given a dataset, methods for statistical causal inference evaluate how an intervention on some treatment variable $A$ would affect a response variable $Y$, typically in the presence of other potential confounding variables. CausalTables.jl is a Julia package for storing and simulating data for causal inference, with three main goals in mind:

1. Providing a `CausalTable` data structure that wraps any dataset stored in a [Tables.jl](https://tables.juliadata.org/stable/)-compatible format with a causal structure (i.e. labeling treatment, response, confounders, and other variables) for use by external causal inference packages.
2. Easily generating random data from a structural causal model (SCM), to be used in simulations.
3. Computing true conditional distributions and causal effect estimands for a given model, allowing users to evaluate the performance of new and existing causal inference methods against the ground truth.



## Installation
CausalTables.jl can be installed using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```
Pkg> add CausalTables
```

## Running a simulation with CausalTables.jl
To simulate data, one must first define a `StructuralCausalModel` (SCM). An SCM is composed of a `DataGeneratingProcess`, which is a sequence of random variables, along with labels for treatment, response, and confounder variables. For example, the following code defines an SCM with a binary treatment $A$, a continuous confounder $W$, and a continuous response $Y$. The `@dgp` macro constructs a `DataGeneratingProcess` object according to the simple syntax `name ~ distribution`, where `rhs` is a `Distribution` object from [Distributions.jl](https://juliastats.org/Distributions.jl/stable/). More advanced syntax is detailed in the [documentation](https://salbalkus.github.io/CausalTables.jl/dev/).

```
using CausalTables
using Distributions

dgp = @dgp(
    W ~ Beta(2, 4),
    A ~ (@. Bernoulli(W)),
    Y ~ (@. Normal(A + W, 1))
)

scm = StructuralCausalModel(dgp; treatment = :A, response = :Y, confounders = [:W])
```

Once a `StructuralCausalModel` is defined, one can then draw a randomly-generated `CausalTable` according to the SCM using the `rand` function:

```
ctbl = rand(scm, 100)
```

Alternatively, it is also possible to approximate the "ground truth" value of a variety of relevant causal estimands from this SCM, including counterfactual means (`cfmean`), as well as average treatment effects (`ate`) and average policy effects (`ape`). For example, the ground truth average treatment effect for this SCM can be approximated like so:

```
ate(scm)

CausalTable
┌───────────┬───────┬───────────┐
│         W │     A │         Y │
│   Float64 │  Bool │   Float64 │
├───────────┼───────┼───────────┤
│  0.144405 │ false │ -0.819133 │
│  0.213987 │ false │    1.2863 │
│ 0.0340754 │ false │  0.917915 │
│  0.203227 │ false │   1.99271 │
│     ⋮     │   ⋮   │     ⋮     │
│  0.140308 │ false │  0.526877 │
│  0.396759 │ false │ -0.308218 │
│  0.217854 │ false │  0.680559 │
│  0.412048 │ false │  0.477174 │
└───────────┴───────┴───────────┘
                  92 rows omitted
Summaries: NamedTuple()
Arrays: NamedTuple()
```

See the [documentation](https://salbalkus.github.io/CausalTables.jl/dev/) for more information and tutorials. 

## Community Guidelines

If you find a bug, have a feature request, or otherwise experience any issues with this software package, please open an issue on the [issue tracker](https://github.com/salbalkus/CausalTables.jl/issues). If you would like to contribute code to the software yourself, we encourage you to open a [pull request](https://github.com/salbalkus/CausalTables.jl/pulls). We welcome all contributions, including bug fixes, documentation improvements, and new features.




