# Network Summaries

Typically, performing causal inference on a network relies on summarizing the treatment and covariates of each unit's neighbors using some sort of **summary function**. For example, in a study evaluating the effect of electric vehicle adoption on air pollution, one might model the commuting patterns between counties as a network and evaluate the effect of the *sum* the number of electric vehicles commuting into each county. The following documents summary measures available in CausalTables.jl.

## Summarizing a CausalTable

```@docs
summarize
```

## Existing Summary Measures

```@docs
NeighborSum
```

```@docs
NeighborProduct
```

```@docs
NeighborMaximum
```

```@docs
NeighborMinimum
```

```@docs
NeighborMode
```

```@docs
Friends
```

## Defining Your Own Summary Measures

