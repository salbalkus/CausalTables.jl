# Approximating Ground Truth Causal Estimands

In causal inference, we are often interested in the value of some causal estimand, such as the average treatment effect (ATE) or the average policy effect (APE). These estimands are defined in terms of counterfactuals $Y(a)$, the value of an outcome $Y$ under a given treatment regime. For example, the average treatment effect is denoted

$$\mathbb{E}\Big(Y(1) - Y(0)\Big)$$

where $Y(a)$ denotes the outcome $Y$ that would have occurred had the treatment been set to $a$. Similarly, we denote the average policy effect as

$$\mathbb{E}\Big(Y(a^*) - Y(a)\Big)$$

where $a$ is the natural value of treatment under no intervention and $a^*$ is the value of treatment under some policy.

CausalTables.jl provides functions that numerically approximate the values of several common estimands given a ground truth `StructuralCausalModel` object. This can be useful for evaluating the performance of causal inference methods on simulated data. Available estimands include:

- Counterfactual Means (`cfmean`)
- Counterfactual Differences (`cfdiff`)
- Average Treatment Effect (`ate`), including among the Treated (`att`) and Untreated (`atu`)
- Average Policy Effect (`ape`), also known as the causal effect of a Modified Treatment Policy.

In addition, one can simulate counterfactual outcomes directly using the `draw_counterfactual` function. Each of these is documented in detail in the following section. For low-level functions that can be used to approximate more complicated custom ground truth estimands in various settings, see [Computing ground truth conditional distributions](ground-truth.md).

## Estimands API

```@autodocs; canonical=false
Modules = [CausalTables]
Order   = [:type, :function]
Pages = ["estimands.jl"]
```

