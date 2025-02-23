---
affiliations:
- index: 1
  name: Department of Biostatistics, Harvard T.H. Chan School of Public
    Health
authors:
- affiliation: 1
  corresponding: true
  name: Salvador V. Balkus
  orcid: 0000-0003-4695-833X
- affiliation: 1
  name: Nima S. Hejazi
  orcid: 0000-0002-7127-2789
bibliography: paper.bib
date: 22 February 2025
tags:
- Julia
- statistics
- causal inference
- tables
title: "CausalTables.jl: Simulating and storing data for statistical
  causal inference in Julia"
toc-title: Table of contents
---

# Summary

Estimating the strength of causal relationships between treatment and
response variables is an important problem across many scientific
disciplines. `CausalTables.jl` is a package that supports causal
inference in Julia by providing two important functionalities. First, it
implements the `CausalTable`, bundling tabular data with a type of
directed acyclic graph (DAG) encoding features' causes. Users can
intervene on treatments and identify causal-relevant variables like
confounders automatically. Second, the package's `StructuralCausalModel`
interface simplifies simulating data from arbitrary causal structures --
and unlike other packages, users can extract ground truth distributions
conditional on the data generated in previous steps. In this way,
`CausalTables.jl` makes it easier to develop and experimentally evaluate
new statistical causal inference methods in Julia.

# Statement of need

The quantitative science of causal inference has emerged over the past
three decades as a set of formalisms for studying cause-and-effect
relationships between variables from observed data
[@pearl2009causality; @hernan2020causal]. Causal inference techniques
have helped scientists and decision-makers better understand important
phenomena in fields ranging from medicine to economics. New software
tools for causal inference are being developed at a rapid pace, but in
the Julia language, there currently do not exist auxiliary tools
designed to support their development. `CausalTables.jl` aims to provide
such a tool.

Implementing and testing causal inference methods in Julia involves two
main challenges. First, causal estimation requires identifying and
modifying features based on their relationships with treatment and
response variables, which might include confounders, mediators, or
instruments. Their required format may differ depending on downstream
packages; for instance, `MLJ.jl` [@blaom2020mlj] requires Table input,
while `GLM.jl` [@bates2023glm] needs a Matrix or `Formula`. Second, when
evaluating a causal estimator on simulated data from a Structural Causal
Model (SCM) [@pearl2009causality], one often desires access to the true
("oracle") conditional distributions of relevant variables in the SCM,
as well as ground truth values of various causal estimands, in order to
test whether the method works correctly.

`CausalTables.jl` addresses both challenges -- the first via the
`CausalTable` interface, which extends `Tables.jl` [@quinn2024tables]
with causal identification routines, and the second via the
`StructuralCausalModel`, which encodes a causal model as a sequence of
conditional distributions from `Distributions.jl`
[@JSSv098i16; @Distributions.jl-2019], providing random sampling and
ground-truth computation. `CausalTables.jl` integrates seamlessly with
established Julia packages, ensuring ease of use for statisticians and
applied scientists alike.

# Comparison to existing packages

While R and Python include many causal packages [@tlverse; @Chen2020],
Julia has relatively fewer. Recent Julia packages for causal inference
include `TMLE.jl` [@TMLE.jl] and `CausalELM.jl` [@CausalELM.jl]. These
focus on specific estimators, rather than general data processing and
simulation like `CausalTables.jl`. The package `CausalInference.jl`
[@Schauer2024] implements causal graphs and discovery algorithms,
similar to CausalDAG [@squires2018causaldag] or DoWhy [@dowhy] in Python
and daggity [@Textor2017] in R. That said, it is generally incompatible
with the tabular data used in practice and does not support simulations.
The simulation capabilities of `CausalTables.jl` are similar to those of
probabilistic programming packages like `Turing.jl` [@turing] or
`Gen.jl` [@gen]. However, while other packages can *sample* data from
SCMs, only `CausalTables.jl` allows extracting *closed-form
distributions* conditional on data drawn in previous steps of the
process.

# Example 1: Data with causal structure

`CausalTables.jl` supports causal inference problems that involve
estimating the effect of at least one treatment on at least one
response. Using the `CausalTable` constructor, one can wrap an existing
Table with causal structure:

::: {.cell execution_count="1"}
``` {.julia .cell-code}
using CausalTables

# Example data in Tables-compatible format
tbl = (W = [0.2, 0.4, 0.7], 
       A = [false, true, true], 
       Y = [0.8, 1.2, 2.3])

# Wrap data as CausalTable
ct_wrap = CausalTable(tbl; treatment = :A, response = :Y, 
                           causes = (A = [:W], Y = [:W, :A]))
```
:::

Convenience functions perform causal data processing. For example, the
general `parents` function selects only features that cause a given
variable; other functions, like `confounders`, select variables with
more specific causal relationships.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
parents(ct_wrap, :Y)
```

::: {.cell-output .cell-output-display execution_count="1"}
    CausalTable
    ┌─────────┬───────┐
    │       W │     A │
    │ Float64 │  Bool │
    ├─────────┼───────┤
    │     0.2 │ false │
    │     0.4 │  true │
    │     0.7 │  true │
    └─────────┴───────┘
    Summaries: NamedTuple()
    Arrays: NamedTuple()
:::
::::

# Example 2: Simulation with ground truth

An SCM defines causal structure by envisaging a data-generating process
as random draws from a sequence of non-parametric structural equations,
with each draw depending on the draws preceding it. For example:

`\begin{align*}
W &\sim Beta(2, 4) \\
A &\sim Bernoulli(0.5W + 0.2) \\
Y &\sim Normal(A + W, 1)
\end{align*}`{=tex}

This SCM can be implemented in `CausalTables.jl` and randomly sampled by
enumerating the sequence of random variables along with labels of their
causal roles:

::: {.cell execution_count="1"}
``` {.julia .cell-code}
using Distributions

# Define sequence of random variables
dgp = @dgp(
    W ~ Beta(2, 4),
    A ~ Bernoulli.(0.5 .* W .+ 0.2),
    Y ~ Normal.(W .+ A, 1)
)

# Define structural causal model
scm = StructuralCausalModel(dgp; 
  treatment = :A, response = :Y
)

ct = rand(scm, 5) # randomly sample
```
:::

Many causal estimands involve applying some intervention to a treatment.
For instance, computing an ATE compares hypothetical responses had
everyone been treated versus no one treated; one can apply these
interventions on a `CausalTable` using the `intervene` function:

::: {.cell execution_count="1"}
``` {.julia .cell-code}
treated = intervene(ct, treat_all)
untreated = intervene(ct, treat_none)
```
:::

After simulating data, the true ("oracle") distribution can be obtained
using `condensity`. Other functions obtain specific features, such as
`conmean` for the conditional mean. These help evaluate how well a
causal estimator might perform if the true distribution were known; for
example, the code below computes the "true" ATE plug-in estimate:

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
mean(conmean(scm, treated, :Y) .- conmean(scm, untreated, :Y))
```

::: {.cell-output .cell-output-display execution_count="1"}
    1.0
:::
::::

`CausalTables.jl` also provides high-level functions to approximate the
ground truth of common causal estimands, such as:

-   Average treatment effects (`ate`) including among the treatment
    (`att`) and untreated (`atu`)
-   Counterfactual means (`cfmean`) and differences (`cfdiff`)
-   Average policy effects (`ape`)

# Closing remarks

The goal of `CausalTables.jl` is to simplify causal inference in Julia.
So far, it has been used to experimentally evaluate novel causal
estimators for continuous treatments on network data [@Balkus2024], and
also been integrated into `TMLE.jl` [@TMLE.jl]. As interest in causal
inference grows, `CausalTables.jl` aims to provide a user-friendly
foundation for practitioners to develop and test new causal methods in
the Julia ecosystem.

# Acknowledgements

Salvador Balkus acknowledges support from the National Institute of
Environmental Health Sciences (award no. T32 ES007142) and the National
Science Foundation (award no. DGE 2140743).

# References {#references .unnumbered}
