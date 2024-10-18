# CausalTables.jl

[![Build Status](https://github.com/salbalkus/CausalTables.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/salbalkus/CausalTables.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/salbalkus/CausalTables.jl/badge.svg?branch=main)](https://coveralls.io/github/salbalkus/CausalTables.jl?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

*A package for storing and simulating data for causal inference in Julia.*

Given a dataset, methods for statistical causal inference evaluate how an intervention on some treatment variable $A$ would affect a response variable $Y$, typically in the presence of other potential confounding variables. CausalTables.jl is a tool for to evaluating and compare the statistical performance of causal inference methods in Julia. The package provides:

1. A `CausalTable` data structure that wraps any dataset stored in a [Tables.jl](https://tables.juliadata.org/stable/)-compatible format with a causal structure (i.e. labeling treatment, response, confounders, and other variables) for use by external causal inference packages.
2. A `StructuralCausalModel` data structure allowing users to easily generate random datasets endowed with a given causal structure.
3. Methods to computing true conditional distributions and approximate ground truth causal effect estimands for a given `StructuralCausalModel`

When used together, both functionalities allow users to benchmark methods for statistical causal inference in Julia -- that is, evaluate the performance of new and existing causal inference methods against both each other and the the ground truth.

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

CausalTable
┌───────────┬───────┬───────────┐
│         W │     A │         Y │
│   Float64 │  Bool │   Float64 │
├───────────┼───────┼───────────┤
│  0.381527 │  true │  0.809227 │
│  0.576206 │  true │   3.22163 │
│  0.380546 │ false │   1.70505 │
│  0.226648 │ false │  0.185022 │
│     ⋮     │   ⋮   │     ⋮     │
│  0.385836 │ false │ -0.392848 │
│  0.204554 │ false │  0.638084 │
│  0.232177 │  true │  0.832707 │
│ 0.0465189 │ false │   1.29168 │
└───────────┴───────┴───────────┘
                  92 rows omitted
Summaries: NamedTuple()
Arrays: NamedTuple()
```

Alternatively, it is also possible to approximate the "ground truth" value of a variety of relevant causal estimands from this SCM, including counterfactual means (`cfmean`), as well as average treatment effects (`ate`) and average policy effects (`ape`). For example, the ground truth average treatment effect for this SCM can be approximated like so:

```
ate(scm)

(μ = 1.0006736394005957, eff_bound = 2.0019300075151616)
```

See the [documentation](https://salbalkus.github.io/CausalTables.jl/dev/) for more information and tutorials. 

## Community Guidelines

If you find a bug, have a feature request, or otherwise experience any issues with this software package, please open an issue on the [issue tracker](https://github.com/salbalkus/CausalTables.jl/issues). If you would like to contribute code to the software yourself, we encourage you to open a [pull request](https://github.com/salbalkus/CausalTables.jl/pulls). We welcome all contributions, including bug fixes, documentation improvements, and new features.




