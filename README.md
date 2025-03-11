# CausalTables.jl

[![Build Status](https://github.com/salbalkus/CausalTables.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/salbalkus/CausalTables.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/salbalkus/CausalTables.jl/badge.svg?branch=main)](https://coveralls.io/github/salbalkus/CausalTables.jl?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![JOSS Status](https://joss.theoj.org/papers/68c43e832d063050a4e67528191e8148/status.svg)](https://joss.theoj.org/papers/68c43e832d063050a4e67528191e8148)

*A common interface for processing and simulating data for causal inference in Julia.*

[Causal inference](https://en.wikipedia.org/wiki/Causal_inference) is the process of estimating, from data, the effect of a treatment variable on an outcome variable -- typically in the presence of confounders. The goal of `CausalTables.jl` is to simplify the development of statistical causal inference methods in Julia. To this end, the package provides two sets of tools:

1. The `CausalTable`, a [Tables.jl](https://tables.juliadata.org/stable/)-compliant data structure that wraps a table of data with labels of the causes of relevant variables, denoted via a type of [directed acyclic graph (DAG)](https://en.wikipedia.org/wiki/Directed_acyclic_graph). Users can call existing functions to easily intervene on treatment variables, identify common subsets of variables (confounders, mediators, instruments, etc.) or use causal labels in other ways -- all while still allowing the data to be used with other Julia packages that accept Tables.jl data structures.
2. The `StructuralCausalModel` interface, which allows users to encode a Structural Causal Model (SCM), a sequence of conditional distributions where each distribution can depend (causally) on any of the previous. This supports simulating data from arbitrary causal structures, extract ground truth distributions conditional on the data generated in previous steps, and approximating common ground-truth estimands such as the average treatment effect or policy effect. 

**What sets this package apart?** `CausalTables.jl` provides a common interface for manipulating tabular data for causal inference. While packages like [CausalInference.jl](https://mschauer.github.io/CausalInference.jl/latest/) only focus on causal graphs and discovery algorithms, the `CausalTable` interface provides utility functions to clean and manipulate practical datasets for input into statistical estimators. The simulation capabilities of `CausalTables.jl` are similar to those of probabilistic programming languages like [Turing.jl](https://turing.ml/dev/) or [Gen.jl](https://www.gen.dev/); however, unlike these packages, with `CausalTables.jl` users can extract the true conditional distributions of relevant variables from a dataset in closed-form *after* data has been generated. This makes it easy to extract parameters like ground-truth ("oracle") conditional means or propensity scores, which are often helpful for testing whether an estimator is behaving as intended.

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

scm = StructuralCausalModel(dgp; treatment = :A, response = :Y, causes = (A = [:W], Y = [:A, :W]))
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

A `CausalTable` is a Table with a causal structure, such as labels for treatment, response, and causes. In addition to implementing the standard Tables.jl interface, CausalTables.jl also provides extra functions to make working with causal data easier. See the [documentation](https://salbalkus.github.io/CausalTables.jl/dev/) for more information.

Given an SCM, it is also possible to approximate the "ground truth" value of a variety of relevant causal estimands from this SCM, including counterfactual means (`cfmean`), as well as average treatment effects (`ate`) and average policy effects (`ape`). For example, the ground truth average treatment effect for this SCM can be approximated like so:

```
ate(scm)

(μ = 1.0, eff_bound = 6.767108201891064)
```

Alternatively, one can compute the ground truth of low-level statistical functionals, such as conditional means or propensity scores, for use in downstream analyses. 

```
propensity(scm, ctbl, :A)

100-element Vector{Float64}:
 0.17900559797871887
 0.653489070854697
 ⋮
 0.8247576289722464
```

See the [documentation](https://salbalkus.github.io/CausalTables.jl/dev/) for more information and tutorials. 

## Community Guidelines

If you find a bug, have a feature request, or otherwise experience any issues with this software package, please open an issue on the [issue tracker](https://github.com/salbalkus/CausalTables.jl/issues). If you would like to contribute code to the software yourself, we encourage you to open a [pull request](https://github.com/salbalkus/CausalTables.jl/pulls). We welcome all contributions, including bug fixes, documentation improvements, and new features.




