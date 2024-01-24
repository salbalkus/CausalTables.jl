# Network Summaries

Typically, performing causal inference on a network relies on summarizing the treatment and covariates of each unit's neighbors using some sort of **summary function**. For example, in a study evaluating the effect of electric vehicle adoption on air pollution, one might model the commuting patterns between counties as a network and evaluate the effect of the *sum* the number of electric vehicles commuting into each county. The following documents summary measures available in CausalTables.jl.

## Summarizing a CausalTable

```@docs
summarize(x::CausalTable, keep = true)
```

## Existing Summary Measures

```@docs
NeighborSum(var_to_summarize::Symbol; use_inneighbors::Bool = true)
```

```@docs
NeighborProduct(var_to_summarize::Symbol; use_inneighbors::Bool = true)
```

```@docs
NeighborMaximum(var_to_summarize::Symbol; use_inneighbors::Bool = true)
```

```@docs
NeighborMinimum(var_to_summarize::Symbol; use_inneighbors::Bool = true)
```

```@docs
NeighborMode(var_to_summarize::Symbol; use_inneighbors::Bool = true)
```

```@docs
Friends(; use_inneighbors::Bool = true)
```

## Defining Your Own Summary Measures

