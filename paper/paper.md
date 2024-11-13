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
date: 2024-11-13
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

Estimating the strength of causal relationships between variables is a
problem of prime importance across scientific disciplines -- and one for
which many competing statistical methods are being developed.
CausalTables.jl provides tools to evaluate and compare statistical
causal inference methods in Julia. The package provides two main
functionalities. First, it implements a `CausalTable` interface for
storing data with partially-labeled causal structure in a
Tables.jl-compatible format. Second it provides a
`StructuralCausalModel` for randomly generating data with a given causal
structure, as well as computing ground truth parameters. When used
together, both functionalities allow users to more easily use and
benchmark the growing number of methods for causal inference in Julia.

# Statement of need

The field of causal inference helps scientists and decision-makers
understand cause-and-effect relationships between variables in data
[@hernan2020causal]. As interest in this field has grown across
disciplines, so too has the development of software tools for estimating
causal effects. While Julia packages for causal inference have begun to
emerge -- including TMLE.jl [@TMLE.jl] and CausalELM.jl [@CausalELM.jl]
for estimation and CausalInference.jl for graph discovery [@Schauer2024]
-- the ecosystem is still in its infancy. Because new methods for causal
inference in various settings are being developed at a rapid pace, it is
important to have tools that make it easy to evaluate and compare their
performance. The goal of CausalTables.jl is to provide such a tool in
Julia.

Currently, those attempting to benchmark causal inference methods in
Julia face two major challenges. First, packages often have inconsistent
interfaces. For example, some packages might require the user to provide
treatment and response variables as individual vectors, while others
might require the entire dataset in a Tables.jl-compliant format, with
treatment and response labeled via strings or symbols. CausalTables.jl
provides a `CausalTable` interface that simplifies packaging data and
auxiliary causal information together.

The second major challenge is that benchmarking often requires
simulating data from a known Structural Causal Model (SCM)
[@pearl2009causality] and comparing estimated effects to some ground
truth. An SCM is a statistical model, typically defined as a sequence of
random draws, with each draw depending on the previous ones. Both the
potential outcomes and graph-based frameworks of causality can be
represented using SCMs. CausalTables.jl provides a simple way for users
to define their own SCM, draw random data from it, and compute or
approximate the true values of several common causal effect parameters.

By addressing these two major challenges, CausalTables.jl helps simplify
and accelerate the development of tools for statistical causal inference
in Julia. The `CausalTable` interface extends Tables.jl, the most common
interface for accessing tabular data in Julia [@quinn2024tables]. The
SCM framework works in conjunction with Distributions.jl, the premier
Julia package for working with random variables
[@JSSv098i16; @Distributions.jl-2019]. By integrating seamlessly with
other common packages in the Julia ecosystem, CausalTables.jl ensures
both compatibility and ease of use for statisticians and students alike.

# Instructional use cases

The standard causal inference problem is to estimate the effect of a
treatment variable $A$ on a response variable $Y$ in the presence of
confounders $W$. One can benchmark causal inference methods in two ways:
either by imposing a causal structure on an existing dataset, or by
drawing new data randomly from a programmatically-defined SCM.

Wrapping an existing dataset with causal structure is easy. The
`CausalTable` constructor creates a Tables-compliant data structure
coupled with causal information about its data. Calling convenience
functions on this object allows users to perform common causal data
processing tasks. For instance, the `responseparents` function can be
used to select only variables upstream from the response.

:::: {.cell execution_count="1"}
::: {.cell-output .cell-output-display execution_count="1"}
    TaskLocalRNG()
:::
::::

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
using CausalTables

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

Simulating causal data for different settings is slightly more involved.
In the remainder of this section, we will present two examples use cases
of how CausalTables.jl can be used as a benchmarking tool.

## Example 1: Average Treatment Effect

The prototypical causal inference problem involves estimating the
average treatment effect (ATE) of a binary treatment $A$ on some outcome
$Y$, in the presence of potential confounders $W$. The ATE describes the
difference in the counterfactual mean of $Y$ had everyone been treated
versus no one treated. An example SCM describing this scenario might be
the following:

`\begin{align*}
W &\sim Beta(2, 4) \\
A &\sim Bernoulli(W) \\
Y &\sim Normal(A + W, 1)
\end{align*}`{=tex}

To compute the ground truth ATE, we can define the SCM above in
CausalTables.jl by defining the sequence of random variables to be drawn
using the `@dgp` macro. Then, we create a `StructuralCausalModel` object
which labels the steps we want to consider as treatment, response, and
confounders.

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
```
:::

We can then randomly draw datasets from the SCM using the `rand`
function, which returns a `CausalTable` object.

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
    │   0.651 │  true │   2.061 │
    │   0.139 │ false │  -0.560 │
    │   0.652 │  true │   2.758 │
    └─────────┴───────┴─────────┘
    Summaries: NamedTuple()
    Arrays: NamedTuple()
:::
::::

At a high level, CausalTables.jl provides functions to approximate
ground truth values of common causal estimands such as the ATE (using
`ate`), along with their efficiency bound -- the lowest possible
variance achievable by a causal estimator of the given quantity.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
ate(scm)
```

::: {.cell-output .cell-output-display execution_count="1"}
    (μ = 1.001, eff_bound = 2.000)
:::
::::

CausalTables.jl also provides a low-level interface allowing users to
(1) apply common interventions to the treatment variable in a
`CausalTable`, and (2) compute ground-truth conditional densities and
functions of them (mean, variance, et cetera) typically used to
construct causal estimators. For example, the code below computes the
difference in conditional means of $Y$ under treatment versus no
treatment.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
ct = rand(scm, 500)
ct_treated = intervene(ct, treat_all)
ct_untreated = intervene(ct, treat_none)

individual_effect = conmean(scm, ct_treated, :Y) .- conmean(scm, ct_untreated, :Y)
mean(individual_effect)
```

::: {.cell-output .cell-output-display execution_count="1"}
    1.000
:::
::::

The above represents the ground-truth plug-in estimate of the individual
treatment effect (outcome regression) for each unit in the dataset.
Alternatively, one can also compute an inverse-probability weighted
(IPW) estimate with ground-truth weights using the `propensity`
function:

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
y = responsematrix(ct) # get the response
a = treatmentmatrix(ct) # get the treatment
ipw = mean(y .* (2 * a .- 1) ./ propensity(scm, ct, :A))
```

::: {.cell-output .cell-output-display execution_count="1"}
    1.160
:::
::::

Finally, as an alternative, one can randomly generate a new
counterfactual response value for each unit in a CausalTable under a
given intervention using `draw_counterfactual`, and compute the ATE
directly:

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
y_treated = draw_counterfactual(scm, ct, treat_all)
y_untreated = draw_counterfactual(scm, ct, treat_none)

mean(y_treated .- y_untreated)
```

::: {.cell-output .cell-output-display execution_count="1"}
    0.951
:::
::::

## Example 2: Modified Treatment Policies

CausalTables.jl can be used for more than just binary treatments; it
also supports more exotic estimands. Consider the following SCM, in
which the treatment $A$ takes on a continuous value.

::: {.cell execution_count="1"}
``` {.julia .cell-code}
dgp = @dgp(
    W1 ~ Poisson(10),
    W2 ~ Bernoulli(0.5),
    A ~ (@. Normal(W1 + 5*W2, 1)),
    Y ~ (@. Normal(2*A + W2*A + 0.5*W1, 1))
)

scm = StructuralCausalModel(dgp; 
  treatment = :A, response = :Y, confounders = [:W1, :W2]
)
```
:::

In the continuous treatment setting, we are often interested in the
effect of a *modified treatment policy* (MTP), which poses the question:
"how would the counterfactual outcome change had we applied some
intervention $d$ to the existing treatment?" [@Haneuse2013]. In this
setting, rather than estimating an ATE, we estimate an *average policy
effect* (APE) -- the difference between $Y$ under the natural treatment
versus the treatment upon which we have intervened. A common example is
the additive MTP $d(a) = a + 1$; when the relationship between $A$ and
$Y$ is linear, this is equivalent to fitting a linear regression, but
using CausalTables.jl, we can obtain approximate the ground truth even
when the relationship is nonlinear.

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
ape(scm, additive_mtp(1))
```

::: {.cell-output .cell-output-display execution_count="1"}
    (μ = 2.499, eff_bound = 2.255)
:::
::::

One strategy for estimating an APE is using a parametric outcome
regression (such as a linear model) that predicts $Y$ under the modified
treatment policy $A = d(a)$ and computes its average difference from the
naturally observed $Y$. We can use CausalTables.jl to see how well this
procedure would work if we knew the true value of the outcome regression
using the `intervene` and `conmean` functions:

:::: {.cell execution_count="1"}
``` {.julia .cell-code}
# Randomly generate data
ct = rand(scm, 1000) 
ct_intervened = intervene(ct, additive_mtp(1))
y = responsematrix(ct)

plugin = mean(conmean(scm, ct_intervened, :Y) .- y)
```

::: {.cell-output .cell-output-display execution_count="1"}
    2.549
:::
::::

# Closing Remarks

The flexibility of CausalTables.jl allows users to easily extract
ground-truth values for any relevant aspect of a data generating
process. This supports benchmarking causal estimators of virtually any
estimand that fits in the SCM framework, not just those mentioned in
these two examples. While the package includes high-level functions to
approximate several common causal estimands, users can also write their
own interventions and use low-level functions such as `intervene`,
`draw_counterfactual`, and `condensity` to approximate the ground truth
of custom causal estimands as well. CausalTables.jl will serve as a
useful tool for statisticians and other scientists seeking to evaluate
various causal inference methods using simulation studies.

# Acknowledgements

Salvador Balkus acknowledges support from the National Institute of
Environmental Health Sciences (award no.\~T32 ES007142) and the National
Science Foundation (award no.\~DGE 2140743).

# References {#references .unnumbered}
