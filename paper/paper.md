---
affiliations:
- index: 1
  name: Harvard T.H. Chan School of Public Health, U.S.A.
authors:
- affiliation: 1
  corresponding: true
  name: Salvador V. Balkus
  orcid: 0000-0003-4695-833X
- affiliation: 1
  name: Nima S. Hejazi
  orcid: 0000-0002-7127-2789
bibliography: paper.bib
date: 2024-11-12
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

CausalTables.jl provides tools to evaluate and compare the statistical
performance of causal inference methods in Julia. The package implements
two main functionalities. First, it defines a `CausalTable` interface
for storing data with partially-labeled causal structure in a
Tables.jl-compatible format. Second, it allows users to define a
structural causal model (SCM) for randomly generating data with a given
causal structure, as well as computing or approximating ground truth
causal effect parameters from that SCM. When used together, both
functionalities allow users to benchmark the growing number of methods
for causal inference in Julia.

# Statement of need

The field of causal inference helps scientists and decision-makers
understand cause-and-effect relationships between variables in data
[@hernan2020causal]. As interest in this field has grown across
disciplines, so too has the development of software tools for estimating
causal effects. In Julia, packages for causal inference have begun to
emerge, such as TMLE.jl [@TMLE.jl] and CausalELM.jl [@CausalELM.jl],
though such packages are generally still in their infancy. Because new
methods for causal inference in various settings are being developed at
a rapid pace, it is important to have tools that make it easy to
evaluate and compare their performance. The goal of CausalTables.jl is
to provide such a tool in Julia.

Currently, those attempting to benchmark causal inference methods in
Julia face two major challenges. First, packages often have inconsistent
interfaces. The canonical problem in causal inference typically takes
the same form across applications: estimate the effect of some treatment
variable $A$ on a response variable $Y$ in the presence of confounders
$W$. Howevever, software packages to do this often require data and
their \`\`causal labels'' to be provided as input in different ways. For
example, some methods might require the user to provide vectors for
treatment and response, while others might require the entire dataset in
a Tables.jl format with treatment and response labels as strings or
symbols. By providing a common interface for storing causal structure
information in a Tables-compatible format, CausalTables.jl makes it easy
to package data and auxiliary causal information and extract the
necessary components needed for benchmarking.

The second major challenge is that evaluating the performance of causal
inference methods typically requires simulating data from a known
Structural Causal Model (SCM) [@pearl2009causality] and comparing
estimated effects to some ground truth value. An SCM is a statistical
model, typically defined as a sequence of draws from probability
distribution. Both potential outcomes and graph-based philosophies of
causality can be represented as an SCM. CausalTables.jl provides a
simple way for users to define their own SCM, draw random datasets from
it, and compute or approximate the true values of several common causal
effect parameters for the SCM at hand.

By providing interfaces to address these two major challenges,
CausalTables.jl will help simplify and accelerate the development of
tools for statistical causal inference in Julia. The `CausalTable`
interface extends Tables.jl, the most common interface for accessing
tabular data in Julia [@quinn2024tables]. The SCM framework works in
conjunction with Distributions.jl, the most popular Julia package for
working with random variables [@JSSv098i16; @Distributions.jl-2019].
Hence, CausalTables.jl integrates seamlessly with other common packages
in the Julia ecosystem, ensuring both compatibility and ease of use for
statisticians and students working in Julia.

# Simulating data with a causal structure

In causal inference, data is conceptualized as being drawn from a
*structural causal model* (SCM). An SCM is typically formulated as a
sequence of draws random variables, with each draw potentially
independent on previous draws. The standard causal inference problem is
to estimate the effect of a treatment variable $A$ on a response
variable $Y$ in the presence of confounders $W$. An example of one such
SCM might be the following:

$$
\begin{align*}
W &\sim Beta(2, 4) \\
A &\sim Bernoulli(W) \\
Y &\sim Normal(A + W, 1)
\end{align*}
$$

We can define the SCM above in CausalTables.jl as follows. First, we
define the sequence of random variables to be drawn using the `@dgp`
macro, which creates a `DataGeneratingProcess`, or DGP. Then, we create
a `StructuralCausalModel` object from the DGP by labeling which steps we
want to consider as treatment, response, and confounders.

::: {.cell execution_count="1"}
``` {.julia .cell-code}
using Distributions
using CausalTables

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
```
:::

We can then randomly draw a dataset from the SCM using the `rand`
function. The function returns a `CausalTable` object, a data structure
implementing the Tables.jl [@quinn2024tables] interface that stores both
the randomly-drawn data, causal structure labels, and additional
metadata output from the DGP, if any.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
ct = rand(scm, 3)
```

::: {.cell-output .cell-output-display execution_count="1"}
    CausalTable
    ┌─────────┬───────┬─────────┐
    │       W │     A │       Y │
    │ Float64 │  Bool │ Float64 │
    ├─────────┼───────┼─────────┤
    │   0.050 │  true │  -0.121 │
    │   0.519 │ false │  -0.002 │
    │   0.301 │ false │   0.356 │
    └─────────┴───────┴─────────┘
    Summaries: NamedTuple()
    Arrays: NamedTuple()
:::
::::

Existing data can also be wrapped as a `CausalTable` using its
constructor. This allows users to store data with partially-labeled
causal structure in a Tables-compatible format as input to external
packages, or use convenience functions to easily partition the table for
causal inference tasks. For instance, the `responseparents` function can
be used to select only variables upstream from the response. An example
is shown below.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
# Example data in a Tables-compatible format
tbl = (W = [0.2, 0.4, 0.7], 
       A = [false, true, true], 
       Y = [0.8, 1.2, 2.3])

# Wrap the data as a CausalTable
ct_wrap = CausalTable(tbl;
                 treatment = :A, 
                 response = :Y, 
                 confounders = [:W])
          
# Select only variables upstream from the response
responseparents(ct_wrap)
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

# Computing ground truth parameters

Causal inference typically involves estimating various parameters under
some intervention on the treatment variable in the original data.
CausalTables.jl provides both high-level and low-level interfaces for
obtaining the ground truth of these parameters.

The low-level interface includes functions to easily (1) apply common
interventions to a CausalTable object, and (2) compute ground-truth
conditional densities and functions of densities (mean, variance, et
cetera).

## Example 1: Average Treatment Effects

The example below computes the difference in conditional means of `Y`
for each unit in the dataset under two interventions: treating all units
and treating none.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
ct = rand(scm, 1000)
ct_treated = intervene(ct, treat_all)
ct_untreated = intervene(ct, treat_none)

individual_effect = conmean(scm, ct_treated, :Y) .- conmean(scm, ct_untreated, :Y)
plugin = mean(individual_effect)
```

::: {.cell-output .cell-output-display execution_count="1"}
    1.000
:::
::::

The above represents the ground-truth plug-in estimate of the individual
treatment effect (outcome regression) for each unit in the dataset.
Alternatively, one can also compute an inverse-probability weighted
estimate with ground-truth weights using the `propensity` function, as
shown in the second example:

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
y = responsematrix(ct) # get the response variable
ipw = mean(y ./ propensity(scm, ct, :A)) .- mean(y)
```

::: {.cell-output .cell-output-display execution_count="1"}
    0.949
:::
::::

Finally, one can randomly generate a new counterfactual response value
for each unit in a CausalTable under a given intervention using
`draw_counterfactual`. The average treatment effect (ATE) is simply the
difference in means of the counterfactual responses under `treat_all`
versus `treat_none` interventions. Analogous treatment effects for
continuous exposures can be computed similarly (for instance, using
`additive_mtp` or `multiplicative_mtp` functions, or by writing custom
intervention functions; see the documentation for more information).

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
y_treated = draw_counterfactual(scm, ct, treat_all)
y_untreated = draw_counterfactual(scm, ct, treat_none)

mean(y_treated .- y_untreated)
```

::: {.cell-output .cell-output-display execution_count="1"}
    1.025
:::
::::

CausalTables.jl also includes a second, high-level interface allowing
users to approximate the ground truth directly for common causal
estimands, such as counterfactual means (`cfmean`), average treatment
effects (`ate`), including among the treated (`att`) and the untreated
(`atu`), and average policy effects (`ape`) under arbitrary
interventions on the treatment. These estimands are approximated using
Monte Carlo integration.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
ate(scm)
```

::: {.cell-output .cell-output-display execution_count="1"}
    (μ = 0.997, eff_bound = 2.001)
:::
::::

## Example 2: Modified Treatment Policies

# Acknowledgements

Salvador Balkus acknowledges support from the National Institute of
Environmental Health Sciences (award no.\~T32 ES007142) and the National
Science Foundation (award no.\~DGE 2140743).

# References {#references .unnumbered}
