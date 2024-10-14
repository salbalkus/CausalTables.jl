# CausalTables.jl

[![Build Status](https://github.com/salbalkus/CausalTables.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/salbalkus/CausalTables.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/salbalkus/CausalTables.jl/badge.svg?branch=main)](https://coveralls.io/github/salbalkus/CausalTables.jl?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Given a dataset, methods for statistical causal inference evaluate how an intervention on some treatment variable $A$ would affect a response variable $Y$, typically in the presence of other potential confounding variables. CausalTables.jl is a Julia package for storing and simulating data for causal inference, with three main goals in mind:

1. Providing a `CausalTable` data structure that wraps any dataset stored in a [Tables.jl](https://tables.juliadata.org/stable/)-compatible format with a causal structure (i.e. labeling treatment, response, confounders, and other variables) for use by external causal inference packages.
2. Easily generating random data from a structural causal model (SCM), to be used in simulations.
3. Computing true conditional distributions and causal effect estimands for a given model, allowing users to evaluate the performance of new and existing causal inference methods against the ground truth.

See the [documentation](https://salbalkus.github.io/CausalTables.jl/dev/) for more information and tutorials. 

## Installation
CausalTables.jl can be installed using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```
Pkg> add CausalTables
```

## Community Guidelines

If you find a bug, have a feature request, or otherwise experience any issues with this software package, please open an issue on the [issue tracker](https://github.com/salbalkus/CausalTables.jl/issues). If you would like to contribute code to the software yourself, we encourage you to open a [pull request](https://github.com/salbalkus/CausalTables.jl/pulls). We welcome all contributions, including bug fixes, documentation improvements, and new features.




