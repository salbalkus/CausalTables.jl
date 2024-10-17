---
title: 'CausalTables.jl: Simulating and storing data for statistical causal inference in Julia'
tags:
  - Julia
  - statistics
  - causal inference
  - tables
authors:
  - name: Salvador V. Balkus
    orcid: 0000-0003-4695-833X
    affiliation: 1
    corresponding: true 
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: 1
affiliations:
 - name: Harvard T.H. Chan School of Public Health, U.S.A.
   index: 1
date: 17 October 2024
bibliography: paper.bib
---

# Summary

CausalTables.jl provides tools to evaluate and compare the statistical performance of causal inference methods in Julia. The package implements two main functionalities. First, it defines a `CausalTable` interface for storing data with partially-labeled causal structure in a Tables.jl-compatible format. Second, it allows users to define a structural causal model (SCM) for randomly generating data with a given causal structure, as well as computing or approximating ground truth causal effect parameters from that SCM. When used together, both functionalities allow users to benchmark the growing number of methods for causal inference in Julia. 

# Statement of need

The field of causal inference helps scientists and decision-makers understand cause-and-effect relationships between variables in data [@hernan2020causal]. As interest in this field has grown across disciplines, so too has the development of software tools for estimating causal effects. In Julia, packages for causal inference have begun to emerge [@TMLE.jl; @CausalELM.jl], though they are generally still in their infancy. Because new methods for causal inference in various settings are being developed at a rapid pace, it is important to have tools that make it easy to evaluate and compare their performance. The goal of CausalTables.jl is to provide such a tool in Julia. 

Currently, those attempting to benchmark causal inference methods in Julia face two major challenges. First, packages often have inconsistent interfaces. The canonical problem in causal inference typically takes the same form across applications: estimate the effect of some treatment variable $A$ on a response variable $Y$ in the presence of confounders $W$. Howevever, software packages to do this often require data and their ``causal labels'' to be provided as input in different ways. For example, some methods might require the user to provide vectors for treatment and response, while others might require the entire dataset in a Tables.jl format with treatment and response labels as strings or symbols. By providing a common interface for storing causal structure information in a Tables-compatible format, CausalTables.jl makes it easy to package data and auxiliary causal information and extract the necessary components needed for benchmarking. 

The second major challenge is that evaluating the performance of causal inference methods typically requires simulating data from a known Structural Causal Model (SCM) [@pearl2009causality] and comparing estimated effects to some ground truth value. An SCM is a statistical model, typically defined as a sequence of draws from probability distribution. Both potential outcomes and graph-based philosophies of causality can be represented as an SCM. CausalTables.jl provides a simple way for users to define their own SCM, draw random datasets from it, and compute or approximate the true values of several common causal effect parameters for the SCM at hand.

By providing interfaces to address these two major challenges, CausalTables.jl will help simplify and accelerate the development of tools for statistical causal inference in Julia. The `CausalTable` interface extends Tables.jl, the most common interface for accessing tabular data in Julia [@quinn2024tables]. The SCM framework works in conjunction with Distributions.jl, the most popular Julia package for working with random variables [@JSSv098i16; @Distributions.jl-2019]. Hence, CausalTables.jl integrates seamlessly with other common packages in the Julia ecosystem, ensuring both compatibility and ease of use for statisticians and students working in Julia. 

# Relationship to existing software

# Example use case

# Acknowledgements

Salvador Balkus acknowledges support from the National Institute of Environmental Health Sciences (award no.~T32 ES007142) and the National Science Foundation (award no.~DGE 2140743).

# References