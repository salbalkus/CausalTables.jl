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
date: 20 January 2025
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
disciplines. `CausalTables.jl` is a Julia package that helps
statisticians and applied scientists create, manipulate, and simulate
datasets labeled with relevant causal structure. Together, its
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
two major challenges. First, statistical causal inference often requires
extracting features from data based on their relationships with
"treatment" and "response" variables; these might include confounders,
mediators, or instruments. The format of these variables might even
differ depending on downstream analysis package; for instance, `MLJ.jl`
[@blaom2020mlj] requires input to be a Table, but `GLM.jl`
[@bates2023glm] necessitates a Matrix or `Formula` object. Second,
testing the performance of new estimators often requires simulating data
for numerical experiments from a Structural Causal Model (SCM)
[@pearl2009causality] so as to compare them to an underlying ground
truth (encoded via interventions on the SCM).

`CausalTables.jl` provides an interface to solve these two problems,
simplifying the development of packages for statistical causal inference
on tabular data in Julia. It implements a `CausalTable` interface that
extends `Tables.jl` [@quinn2024tables] to also store necessary causal
relationships between variables. In addition, the package implements a
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
estimation techniques, each implementing different interfaces to label
causal structure for their specific problems. They do not provide a
general simulation or causal-specific data processing framework like
`CausalTables.jl`. On the other hand, `CausalInference.jl`
[@Schauer2024] provides an interface for representing causal graphs and
implements causal discovery algorithms, similar to CausalDAG
[@squires2018causaldag] or DoWhy [@dowhy] in Python and daggity
[@Textor2017] in R. However, it is generally incompatible with the
tabular data format required by statistical tools, and also cannot
simulate data. In fact, as far as we are aware, `CausalTables.jl` is the
first package for simulating and extracting ground-truth causal
estimands from an existing SCM in Julia.

# Example 1: Data Preprocessing

`CausalTables.jl` supports causal inference problems that involve
estimating the effect of at least one treatment variable $A$ on a
response variable $Y$. Using the `CausalTable` constructor, one can wrap
existing data as a `Tables.jl`-compliant structure coupled with causal
structure labels.

::: {.cell execution_count="1"}
``` {.julia .cell-code}
using CausalTables

# Example data in Tables-compatible format
tbl = (W = [0.2, 0.4, 0.7], 
       A = [false, true, true], 
       Y = [0.8, 1.2, 2.3])

# Wrap data as CausalTable
ct_wrap = CausalTable(tbl; treatment = :A, response = :Y)
```
:::

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
    │     0.2 │ false │
    │     0.4 │  true │
    │     0.7 │  true │
    └─────────┴───────┘
    Summaries: NamedTuple()
    Arrays: NamedTuple()
:::
::::

## Example 2: Simulating data with ground-truth approximations

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
  treatment = :A, response = :Y
)

ct = rand(scm, 500) # randomly draw from the SCM
```
:::

`CausalTables.jl` provides high-level functions that approximate ground
truth values of common causal estimands when called on the `scm`. These
include:

-   Average treatment effects (`ate`) including among the treatment
    (`att`) and untreated (`atu`)
-   Counterfactual means (`cfmean`) and differences (`cfdiff`)
-   Average policy effects (`ape`)

In addition, `CausalTables.jl` implements low-level interface for (1)
applying common interventions to the treatment variable in a
`CausalTable`, (2) drawing randomly from counterfactual distributions,
and (3) computing ground truth conditional densities and functions of
these (e.g., means, variances, propensity scores), which often arise in
the definition of many estimands.. For example, below we compute the
difference in the conditional mean of $Y$ under treatment versus no
treatment, the difference of which is the ATE.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
treated = intervene(ct, treat_all)    # CausalTable with everyone treated
untreated = intervene(ct, treat_none) # CausalTable with no one treated
mean(conmean(scm, treated, :Y) .- conmean(scm, untreated, :Y))
```

::: {.cell-output .cell-output-display execution_count="1"}
    1.0
:::
::::

# Closing remarks

Not only does `CausalTables.jl` provide high-level functions for common
data processing and simulation tasks in causal inference, it can also be
easily extended to support more novel methods and estimands using
low-level functions. The `CausalTable` stores all relevant causal
relationships needed to extract variables related to treatment and
response variables. The `StructuralCausalModel` support simulating data
from any SCM that can be expressed as a sequence of random variables.
Hence, `CausalTables.jl` serves as a useful tool in Julia for both
developing new methods and providing input to existing ones.

# Acknowledgements

Salvador Balkus acknowledges support from the National Institute of
Environmental Health Sciences (award no. T32 ES007142) and the National
Science Foundation (award no. DGE 2140743).

# References {#references .unnumbered}
