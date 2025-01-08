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
date: 08 January 2025
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

Estimating the strength of causal relationships between variables is an
important problem across many scientific disciplines. `CausalTables.jl`
provides tools to (1) easily store and process tabular data endowed with
causal structure and (2) simulate data from causal models for
experimental testing and compute ground-truth estimates. Together, these
functionalities expand the Julia ecosystem by supporting the development
and experimental assessment of new statistical causal inference methods
in Julia.

# Statement of need

The quantitative science of causal inference has emerged over the past
three decades as a set of formalisms for studying cause-and-effect
relationships between variables from observed data
[@pearl2009causality; @hernan2020causal]. Causal inference techniques
have helped applied scientists and decision-makers better understand
important phenomena in fields ranging from health and medicine to
politics and economics. New software tools for causal inference are
being developed at a rapid pace, but in the Julia language, there
currently do not exist auxiliary tools designed to support their
development. `CausalTables.jl` aims to provide such a tool.

Attempts to implement and test causal inference methods in Julia face
two major challenges. First, statistical causal inference requires data
to be transformed in various ways based on the underlying causal
structure, in order to be provided as input to other packages. For
instance, causal methods using `MLJ.jl` [@blaom2020mlj] would
necessitate extracting the variables hypothesized to cause the treatment
and response as Tables; meanwhile; using `GLM.jl` [@bates2023glm] would
require the same but as Matrix or `Formula` objects. Second, testing the
performance of new estimators often requires simulating data for
numerical experiments from a Structural Causal Model (SCM)
[@pearl2009causality] so as to compare them to an underlying ground
truth (encoded via interventions on the SCM).

`CausalTables.jl` provides an interface to solve these two problems,
simplifying the development of packages for statistical causal inference
on tabular data in Julia. It implements a `CausalTable` interface that
extends `Tables.jl`, the most common interface for accessing tabular
data in Julia [@quinn2024tables]. The package also implements a
`StructuralCausalModel` interface for sampling from any SCM and
computing ground-truth estimates of causal parameters. This interface
operates in conjunction with `Distributions.jl`, the primary Julia
package for working with random variables
[@JSSv098i16; @Distributions.jl-2019]. By integrating seamlessly with
other commonly used packages in the Julia ecosystem, `CausalTables.jl`
ensures both compatibility and ease of use for statisticians and applied
scientists alike.

# Comparison to existing packages

While the R and Python ecoysystems include many implementations of
causal methods [@tlverse; @Chen2020], Julia has relatively fewer. Recent
Julia packages for causal inference include `TMLE.jl` [@TMLE.jl] and
`CausalELM.jl` [@CausalELM.jl]. These packages focus on specific
estimation techniques using tabular data, each implementing different
interfaces to label causal structure for their specific causal problems;
they do not provide a general simulation or causal-specific data
processing framework like `CausalTables.jl`. On the other hand,
`CausalInference.jl` [@Schauer2024] provides an interface for
representing causal graphs and implements causal discovery algorithms,
similar to CausalDAG [@squires2018causaldag] or DoWhy [@dowhy] in Python
and daggity [@Textor2017] in R. However, it is generally incompatible
with the tabular data format required by statistical tools, and also
cannot simulate data. In fact, as far as we are aware, `CausalTables.jl`
is the first package for simulating and extracting ground-truth causal
estimands from an existing SCM in Julia.

# Example 1: Data Preprocessing

`CausalTables.jl` supports causal inference problems that involve
estimating the effect of at least one treatment variable $A$ on a
response variable $Y$ in the presence of confounders $W$. Using the
`CausalTable` constructor, one can wrap existing data as a
`Tables.jl`-compliant structure coupled with causal structure labels.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
using CausalTables

# Example data in a Tables-compatible format
tbl = (W = [0.2, 0.4, 0.7], 
       A = [false, true, true], 
       Y = [0.8, 1.2, 2.3])

# Wrap the data as a CausalTable
ct_wrap = CausalTable(tbl; treatment = :A, response = :Y, confounders = [:W])
```

::: {.cell-output .cell-output-display execution_count="1"}
    CausalTable
    ┌─────────┬───────┬─────────┐
    │       W │     A │       Y │
    │ Float64 │  Bool │ Float64 │
    ├─────────┼───────┼─────────┤
    │   0.200 │ false │   0.800 │
    │   0.400 │  true │   1.200 │
    │   0.700 │  true │   2.300 │
    └─────────┴───────┴─────────┘
    Summaries: NamedTuple()
    Arrays: NamedTuple()
:::
::::

Convenience functions perform data processing tasks common to causal
inference, such as selecting or intervening on specific variables. For
example, the `parents` function can be used to select only variables
denoted as causes of $Y$:

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
    │   0.200 │ false │
    │   0.400 │  true │
    │   0.700 │  true │
    └─────────┴───────┘
    Summaries: NamedTuple()
    Arrays: NamedTuple()
:::
::::

## Example 2: Simulating data with ground-truth ATE

An SCM defines causal structure by envisaging a data-generating process
as random draws from a sequence of non-parametric structural equations,
with each draw depending on realizations from draws preceding it. An
example is the following:

`\begin{align*}
W &\sim Beta(2, 4) \\
A &\sim Bernoulli(W) \\
Y &\sim Normal(A + W, 1)
\end{align*}`{=tex}

Using `CausalTables.jl`, we can define the SCM above by enumerating the
sequence of random variables to be drawn, along with labels for the role
they play in estimation, and then randomly draw from it.

::: {.cell execution_count="1"}
``` {.julia .cell-code}
using Distributions

# Define the sequence of random variables to be drawn
dgp = @dgp(
    W ~ Beta(2, 4),
    A ~ Bernoulli.(0.5 .* W .+ 0.2),
    Y ~ Normal.(W .+ A, 1)
)

# Create a structural causal model (SCM) from the DGP
scm = StructuralCausalModel(dgp; 
  treatment = :A, response = :Y, confounders = [:W]
)

ct = rand(scm, 500) # randomly draw from the SCM
```
:::

`CausalTables.jl` provides high-level functions to approximate ground
truth values of common causal estimands, including:

-   Average Treatment Effects (ATE) including among the treatment (ATT)
    and untreated (ATT)
-   Counterfactual Means and Differences
-   Average Policy Effects (APE)

For example, we can compute the ATE on the SCM above like so:

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
ate(scm) # average treatment effect
```

::: {.cell-output .cell-output-display execution_count="1"}
    (μ = 1.000, eff_bound = 2.000)
:::
::::

In addition, `CausalTables.jl` provides a low-level interface allowing
users to (1) apply common interventions to the treatment variable in a
`CausalTable`, (2) draw randomly from counterfactual distributions, and
(3) compute ground truth conditional densities and functions of these
(e.g., mean, variance, propensity scores), which typically arise as
nuisance parameters in the construction of estimators in causal
inference. For example, below, we compute the difference in the
conditional mean of $Y$ under treatment versus no treatment, the
difference of which is the ATE.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
treated = intervene(ct, treat_all)    # CausalTable with everyone treated
untreated = intervene(ct, treat_none) # CausalTable with no one treated
mean(conmean(scm, treated, :Y) .- conmean(scm, untreated, :Y))
```

::: {.cell-output .cell-output-display execution_count="1"}
    1.000
:::
::::

# Closing remarks

`CausalTables.jl` provides useful auxiliary functions to support causal
inference methods on tabular data in Julia that involve one or more
treatment variables and responses. Users can simulate data from any SCM
and benchmark methods using either high-level functions for common
estimands or low-level functions for more exotic estimands. By combining
this with the power of the `CausalTable` interface for processing data
once it is generated, `CausalTables.jl` serves as a useful tool for
scientists seeking to develop and experimentally evaluate new causal
inference methods.

# Acknowledgements

Salvador Balkus acknowledges support from the National Institute of
Environmental Health Sciences (award no. T32 ES007142) and the National
Science Foundation (award no. DGE 2140743).

# References {#references .unnumbered}
