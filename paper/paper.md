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
date: 2025-01-07
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
supports the development of new statistical methods for causal inference
in Julia by providing tools to (1) easily store and process tabular data
endowed with causal structure and (2) simulate data from causal models
for experimental testing. Firstly, the package implements a
`CausalTable` structure that stores features annotated with causal
labels in a `Tables.jl`-compatible format. Its interface includes
causal-relevant functions, such as extracting relevant variables and
applying interventions on treatment. Secondly, `CausalTables.jl`
introduces a `StructuralCausalModel` for randomly generating data from
user-specified causal models and computing ground truth parameters under
the given experiment. Together, these functionalities expand the Julia
ecosystem by supporting the development and experimental assessment of
the growing number of causal inference methods.

# Statement of need

The quantitative science of causal inference has emerged over the past
three decades as a set of formalisms for studying cause-and-effect
relationships between variables from observed data
[@pearl2009causality; @hernan2020causal]; causal inference techniques
have helped applied scientists and decision-makers better understand
important phenomena in fields ranging from health and medicine to
politics and economics. New methods for causal inference are being
developed at a rapid pace, but there currently do not exist auxiliary
tools designed to support their development in the Julia language.
`CausalTables.jl` aims to provide such a tool. Presently, attempts to
implement and test causal inference methods in Julia face two major
challenges.

First, causal inference requires data to be preprocessed in various ways
based on the underlying causal structure. Suppose one were to write
their own method by building on existing statistical packages in Julia.
Using `MLJ.jl` [@blaom2020mlj] would necessitate extracting the
treatment and response as Vectors and the variables hypothesized to
cause them as Tables; meanwhile; using `GLM.jl` [@bates2023glm] would
require the same but as Matrix or `Formula` objects. This challenge has
also led to differences in the API for existing causal methods: for
instance, `CausalELM.jl` [@CausalELM.jl] requires the user to split
apart treatment and response variables as individual vectors, while
`TMLE.jl` require the entire dataset in a `Tables.jl`-compliant format
with treatment and response variables being labeled via strings or
symbols. `CausalTables.jl` provides a `CausalTable` interface that, by
packaging the data and auxiliary causal knowledge together, allows
extracting relevant causal components in multiple ways. This simplifies
both writing new packages as well as processing data to evaluate
existing packages.

Second, testing the performance of new estimators often requires
simulating data for numerical experiments from a Structural Causal Model
(SCM) [@pearl2009causality] so as to compare them to an underlying
ground truth (encoded via interventions on the SCM). An SCM defines
causal structure by envisaging a data-generating process as random draws
from a sequence of non-parametric structural equations, with each draw
depending on realizations from draws preceding it. `CausalTables.jl`
provides a simple, user-friendly way to define an SCM, sample data
randomly from it, and compute or approximate the underlying true values
of several common causal effect parameters.

By addressing these two major challenges---preprocessing and
simulation--- `CausalTables.jl` simplifies and accelerates the
development of tools for statistical causal inference on tabular data in
Julia. The `CausalTable` interface extends `Tables.jl`, the most common
interface for accessing tabular data in Julia [@quinn2024tables]. The
SCM framework operates in conjunction with `Distributions.jl`, the
primary Julia package for working with random variables
[@JSSv098i16; @Distributions.jl-2019]. By integrating seamlessly with
other commonly used packages in the Julia ecosystem, `CausalTables.jl`
ensures both compatibility and ease of use for statisticians and applied
scientists alike.

# Comparison to existing packages

As interest in causal inference continues to grow across disciplines, so
too has the development of software tools for estimating causal effects.
While the a multitude of methods have been implemented in the R and
Python languages (for instance, [@tlverse] or [@Chen2020]), Julia has
seen relatively fewer. Recent Julia packages for causal inference
include `TMLE.jl` [@TMLE.jl] and `CausalELM.jl`[@CausalELM.jl]. These
packages focus on estimation techniques using tabular data: they
implement specific ways to label causal structure, but do not provide a
general simulation or causal-specific data processing interface like
`CausalTables.jl`. On the other hand, `CausalInference.jl`
[@Schauer2024] provides an interface for representing causal graphs and
implements causal discovery algorithms, similar to CausalDAG
[@squires2018causaldag] or DoWhy [@dowhy] in Python and daggity
[@Textor2017] in R. However, it is generally incompatible with the
tabular data format required by statistical tools, and also cannot
simulate data. In fact, as far as we are aware, `CausalTables.jl` is the
first package for simulating and extracting ground-truth causal
estimands from an existing SCM in Julia.

# Instructional use cases

A standard causal inference problem is to estimate the effect of one
treatment variable $A$ on a response variable $Y$ in the presence of
confounders $W$. One can evaluate the performance of causal inference
methods in two ways: either by imposing a causal structure on an
existing dataset, or by drawing new data randomly from a
programmatically-defined SCM.

Wrapping an existing dataset with causal structure is easy. The
`CausalTable` constructor creates a `Tables.jl`-compliant data structure
that wraps any existing data that already satisfies the `Tables.jl`
interface with additional information about its causal structure.
Calling convenience functions on this object allows users to perform
data processing tasks common in causal inference, such as selecting or
intervening on specific variables. For example, the `parents` function
can be used to select only variables denoted as causes of another given
variable.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
using CausalTables

# Example data in a Tables-compatible format
tbl = (W = [0.2, 0.4, 0.7], 
       A = [false, true, true], 
       Y = [0.8, 1.2, 2.3])

# Wrap the data as a CausalTable
ct_wrap = CausalTable(tbl; treatment = :A, response = :Y, confounders = [:W])
          
# Select only variables upstream from the response
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

Simulating causal data for different settings is slightly more involved.
In the remainder of this section, we will present two example use cases
of how `CausalTables.jl` can be used as a benchmarking tool.

## Example 1: Average Treatment Effect

The prototypical causal inference problem involves estimating the
average treatment effect (ATE) of a binary treatment $A$. The ATE
describes the difference in the counterfactual mean of $Y$ had every
unit been treated versus no unit treated. An example SCM describing this
scenario might be the following: `\begin{align*}
W &\sim Beta(2, 4) \\
A &\sim Bernoulli(W) \\
Y &\sim Normal(A + W, 1)
\end{align*}`{=tex} To compute the ground truth ATE via
`CausalTables.jl`, we define the SCM above by enumerating the sequence
of random variables to be drawn using the `@dgp` macro. Then, we create
a `StructuralCausalModel` object which labels the steps we want to
consider as treatment, response, and confounders. Finally, we randomly
draw datasets from the newly instantiated `StructuralCausalModel` using
the `rand` function.

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

`CausalTables.jl` provides functions to approximate ground truth values
of common causal estimands such as the ATE (using `ate`). These
functions also simultaneously estimate the corresponding efficiency
bound---the asymptotic lower bound on the variance---for a class of
estimators (those that are regular and asymptotically linear) commonly
used in causal inference. This facilitates the comparison of candidate
estimators not only in terms of their bias (average distance from the
ground truth) but also their efficiency. Below, we demonstrate these for
the example SCM given above.

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
`CausalTable`, and (2) compute ground truth conditional densities and
functions of these (e.g., mean, variance), which typically arise as
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

The above recovers an estimate of the ground truth via plug-in estimates
based on the outcome regression (i.e., the conditional expectation of
the response, given treatment and covariates). Alternatively, one can
also compute an inverse probability weighted (IPW) estimate with ground
truth weights using the `propensity` function:

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
y = responsematrix(ct) # get the response
a = treatmentmatrix(ct) # get the treatment
mean(y .* (2 * a .- 1) ./ propensity(scm, ct, :A))
```

::: {.cell-output .cell-output-display execution_count="1"}
    0.986
:::
::::

Finally, as an alternative, one can randomly generate a new
counterfactual response value for each unit in a `CausalTable` under a
given intervention using `draw_counterfactual`, and then compute the ATE
directly:

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
y_treated = draw_counterfactual(scm, ct, treat_all)
y_untreated = draw_counterfactual(scm, ct, treat_none)
mean(y_treated .- y_untreated)
```

::: {.cell-output .cell-output-display execution_count="1"}
    1.075
:::
::::

## Example 2: Modified Treatment Policies

`CausalTables.jl` is not limited to settings with binary treatments; it
also supports other estimands. Consider the following SCM, in which the
treatment $A$ is continuous-valued.

::: {.cell execution_count="1"}
``` {.julia .cell-code}
dgp = @dgp(
    W1 ~ Poisson(10),
    W2 ~ Bernoulli(0.5),
    A ~ (@. Normal(W1 + 5*W2, 1)),
    Y ~ (@. Normal(2*A + W2*A + 0.5*W1, 1 + W2))
)

scm = StructuralCausalModel(dgp; 
  treatment = :A, response = :Y, confounders = [:W1, :W2]
)
```
:::

In the continuous treatment setting, a common causal estimand is the
effect of a *modified treatment policy* (MTP), which corresponds to the
question: "how would the counterfactual response change had an
intervention $d(a, w; \delta)$ been applied to the observed treatment
$a$?" [@Haneuse2013]. In this setting, rather than estimating an ATE, we
estimate an *average policy effect* (APE)---the difference between $Y$
under the natural treatment value $a$ and the intervened-upon treatment
that results from the MTP $d(a, w; \delta)$. A common example is the
additive MTP $d(a, w; \delta) = a + \delta$; when the relationship
between $A$ and $Y$ is known to be linear, this is equivalent to the
slope in a linear regression model, but, using `CausalTables.jl`, we can
approximate the ground truth APE even when the relationship is
nonlinear. We demonstrate this below for an MTP indexed by the choice
$\delta = 1$.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
ape(scm, additive_mtp(1)) # average policy effect
```

::: {.cell-output .cell-output-display execution_count="1"}
    (μ = 2.500, eff_bound = 5.252)
:::
::::

One strategy for estimating an APE is to fit a parametric outcome
regression model and use it to predict the outcome $Y$ under the
modified treatment policy $d(A, W;\delta): A \to A + \delta$; the
average difference of these predictions on the observed data yields the
APE. We can use `CausalTables.jl` to simulate the output of such a
procedure had we known the true value of the outcome regression:

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
ct = rand(scm, 500)  # Randomly draw data
ct_intervened = intervene(ct, additive_mtp(1))  # apply MTP
mean(conmean(scm, ct_intervened, :Y) .- responsematrix(ct))
```

::: {.cell-output .cell-output-display execution_count="1"}
    2.528
:::
::::

# Closing remarks

`CausalTables.jl` provides useful auxiliary functions to support causal
inference methods on tabular data in Julia. The package focuses on tools
relevant to estimating the effect of one or more treatment variables on
a response. The `StructuralCausalModel` allows users to easily extract
ground truth values for any relevant aspect of a data-generating
process, supporting the benchmarking of many common causal inference
methods. While the package includes high-level functions to approximate
several prominent estimands, users can also write their own
interventions and use low-level functions such as `intervene`,
`draw_counterfactual`, and `condensity` to approximate the ground truth
of novel estimands. By combining this with the power of the
`CausalTable` interface for processing data once it is generated,
`CausalTables.jl` serves as a useful tool for scientists seeking to
develop and experimentally evaluate new causal inference methods.

# Acknowledgements

Salvador Balkus acknowledges support from the National Institute of
Environmental Health Sciences (award no. T32 ES007142) and the National
Science Foundation (award no. DGE 2140743).

# References {#references .unnumbered}
