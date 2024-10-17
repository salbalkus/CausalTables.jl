# Network Summaries

In network causal inference, methods often rely on summarizing the treatment and covariates of each unit's neighbors using some sort of **summary function**. For example, in a study evaluating the effect of electric vehicle adoption on air pollution, one might model the commuting patterns between counties as a network and evaluate the effect of the *sum* the number of electric vehicles commuting into each county. 

This section documents all available summary measures and how to summarize them within a CausalTable.

## Summarizing a CausalTable

Data wrapped in a CausalTable includes a NamedTuple `summaries` which describes extra variables represented as summary variables over the network. These summary measures can be computed and added to the table by calling the `summarize` function on the CausalTable object.

```@docs; canonical=false
summarize
```

## Existing Summary Measures

The following lists summary measures currently available off-the-shelf in CausalTables.jl. Examples on their use are provided in [Generating Data for Statistical Experiments](generating-data.md) and [Turning data into a `CausalTable`](formatting.md).

```@docs; canonical=false
Sum
Mean
AllOrderStatistics
KOrderStatistics
```


## Defining Your Own Summary Measures

Forthcoming.

