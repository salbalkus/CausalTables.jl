# Approximating Ground Truth Causal Estimands

In causal inference, we are often interested in the value of some causal estimand, such as the average treatment effect (ATE) or the average policy effect (APE). These estimands are defined in terms of counterfactual outcomes, which are the outcomes that would have been observed had the treatment been different. CausalTables.jl provides functions that numerically approximate the values of several common estimands given a ground truth 
`StructuralCausalModel` object. This can be useful for evaluating the performance of causal inference methods on simulated data. 

## Estimands API

```@autodocs
Modules = [CausalTables]
Order   = [:type, :function]
Pages = = ["estimands.jl"]
```

