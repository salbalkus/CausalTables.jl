# Network Summaries

Typically, performing causal inference on a network relies on summarizing the treatment and covariates of each unit's neighbors using some sort of **summary function**. For example, in a study evaluating the effect of electric vehicle adoption on air pollution, one might model the commuting patterns between counties as a network and evaluate the effect of the *sum* the number of electric vehicles commuting into each county. The following documents summary measures available in CausalTables.jl and how to summarize a CausalTable

## Summarizing a CausalTable

Data wrapped in a CausalTable includes a NamedTuple `summaries` which describes extra variables represented as summary variables over the network. These summary measures can be computed and added to the table by calling the `summarize` function on the CausalTable object.

```@docs
summarize
```

## Existing Summary Measures

The following lists summary measures currently available off-the-shelf in CausalTables.jl. Examples on their use are provided in [Generating Data for Statistical Experiments](generating-data.md) and [Turning data into a `CausalTable`](formatting.md).

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

If you want to use a custom summary measure, you can define your own by creating a new type that is a subtype of `NetworkSummary`. The following example shows how to define a new summary measure that computes the mode of a variable over a unit's neighbors:

```@example
mutable struct NeighborMode <: NetworkSummary 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    NeighborMode(var_to_summarize::Symbol; use_inneighbors::Bool = true) = new(var_to_summarize, use_inneighbors)
end
```

Then, simply add a new function to the `summarize` method by defining `summarize(x::CausalTable, summary::NetworkSummary)` function, but replacing the `NetworkSummary` with your new summary type. One easy way to define a new `summarize` function is to use `apply_function_over_neighbors`. This takes a function taking a vector as input and applies it to the variable specified in second argument, filtered to each unit's neighbors. 

In this case, we can use `StatsBase.mode` to compute the mode of each unit's neighbors. 

```@docs
apply_function_over_neighbors
```

```@example
summarize(x::CausalTable, summary::NeighborMode) = apply_function_over_neighbors(x, summary.var_to_summarize, StatsBase.mode; use_inneighbors = summary.use_inneighbors)
```

