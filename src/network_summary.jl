"""
    abstract type NetworkSummary

Abstract type representing a summary of a network.

"""
abstract type NetworkSummary end

get_var_to_summarize(x::NetworkSummary) = x.var_to_summarize


mutable struct NeighborSum <: NetworkSummary 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    
    @doc raw"""
        NeighborSum(var_to_summarize::Symbol; use_inneighbors::Bool = true)

    Constructs a Network Summary denoting the **sum** of `var_to_summarize` over each unit's neighbors. If `var_to_summarize` is ``X``, then mathematically this computes

    ```math
    \sum_{j \in \mathcal{F}_i} X_j
    ```

    where ``F_i`` represents the "friends" or neighbors of unit ``i``.

    # Arguments
    - `var_to_summarize::Symbol`: The variable to summarize.
    - `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

    # Returns
    - `NeighborSum`: The constructed NeighborSum object.
    """
    NeighborSum(var_to_summarize::Symbol; use_inneighbors::Bool = true) = new(var_to_summarize, use_inneighbors)
end

mutable struct NeighborProduct <: NetworkSummary 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    @doc raw"""
        NeighborProduct(var_to_summarize::Symbol; use_inneighbors::Bool = true)

    Constructs a Network Summary denoting the **product** of `var_to_summarize` over each unit's neighbors. If `var_to_summarize` is ``X``, then mathematically this computes

    ```math
    \prod_{j \in \mathcal{F}_i} X_j
    ```

    where ``F_i`` represents the "friends" or neighbors of unit ``i``.

    # Arguments
    - `var_to_summarize::Symbol`: The variable to summarize.
    - `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

    # Returns
    - `NeighborProduct`: The constructed NeighborProduct object.
    """
    NeighborProduct(var_to_summarize::Symbol; use_inneighbors::Bool = true) = new(var_to_summarize, use_inneighbors)
end

abstract type NeighborOrderStatistic <: NetworkSummary end

mutable struct NeighborMaximum <: NeighborOrderStatistic 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    @doc raw"""
        NeighborMaximum(var_to_summarize::Symbol; use_inneighbors::Bool = true)

    Constructs a Network Summary denoting the **maximum** of `var_to_summarize` over each unit's neighbors. If `var_to_summarize` is ``X``, then mathematically this computes

    ```math
    \max_{j \in \mathcal{F}_i} X_j
    ```

    where ``F_i`` represents the "friends" or neighbors of unit ``i``.

    # Arguments
    - `var_to_summarize::Symbol`: The variable to summarize.
    - `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

    # Returns
    - `NeighborMaximum`: The constructed NeighborMaximum object.
    """
    NeighborMaximum(var_to_summarize::Symbol; use_inneighbors::Bool = true) = new(var_to_summarize, use_inneighbors)
end

mutable struct NeighborMinimum <: NeighborOrderStatistic 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    @doc raw"""
        NeighborMinimum(var_to_summarize::Symbol; use_inneighbors::Bool = true)

    Constructs a Network Summary denoting the **minimum** of `var_to_summarize` over each unit's neighbors. If `var_to_summarize` is ``X``, then mathematically this computes

    ```math
    \min_{j \in \mathcal{F}_i} X_j
    ```

    where ``F_i`` represents the "friends" or neighbors of unit ``i``.

    # Arguments
    - `var_to_summarize::Symbol`: The variable to summarize.
    - `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

    # Returns
    - `NeighborMinimum`: The constructed NeighborMinimum object.
    """
    NeighborMinimum(var_to_summarize::Symbol; use_inneighbors::Bool = true) = new(var_to_summarize, use_inneighbors)
end

mutable struct NeighborMode <: NetworkSummary 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    @doc raw"""
        NeighborMode(var_to_summarize::Symbol; use_inneighbors::Bool = true)

    Constructs a Network Summary denoting the **mode** of `var_to_summarize` over each unit's neighbors -- that is, the value occuring most often among units connected to each unit in the network.

    # Arguments
    - `var_to_summarize::Symbol`: The variable to summarize.
    - `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

    # Returns
    - `NeighborMode`: The constructed NeighborMode object.
    """
    NeighborMode(var_to_summarize::Symbol; use_inneighbors::Bool = true) = new(var_to_summarize, use_inneighbors)
end

mutable struct Friends <: NetworkSummary 
    use_inneighbors::Bool
    @doc raw"""
        Friends(var_to_summarize::Symbol; use_inneighbors::Bool = true)

    Constructs a Network Summary denoting the **number of neighbors** of each unit in the network.
    # Arguments
    - `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

    # Returns
    - `Friends`: The constructed Friends object.
    """
    Friends(; use_inneighbors::Bool = true) = new(use_inneighbors)
end


"""
    summarize(x::CausalTable, keep = true)

Summarize a CausalTable by computing the summaries specified in the `summaries` field of the CausalTable and adding them to the existing data.

# Arguments
- `x::CausalTable`: The CausalTable to be summarized.
- `keep_original::Bool`: Whether to return a CausalTable with summarized variables appended (`true`) or a Tables-compatible NamedTuple with only the summarized variables (`false`).

# Returns
- If `keep_original` is `true`, a new CausalTable with merged columns, treatment, response, graph, and summaries. Otherwise, NamedTuple of just variable summaries

# Example
```jldoctest
using Graphs

# Construct a CausalTable
ctbl = CausalTable(
    (A = [1, 2, 3],);
    graph = Graphs.complete_graph(3),
    summaries = (As = NeighborSum(:A), Am = NeighborMaximum(:A))
)

# Compute the summaries of the CausalTable
summarize(ctbl)

# output
CausalTable((A = [1, 2, 3], As = [5, 4, 3], Am = [3, 3, 2]), nothing, nothing, nothing, SimpleGraph{Int64}(3, [[2, 3], [1, 3], [1, 2]]), (As = NeighborSum(:A, true), Am 
= NeighborMaximum(:A, true)))
```
"""
function summarize(x::CausalTable; keep_original = true)
    sumtbl = (;zip(propertynames(x.summaries),  [summarize(x, s) for s in x.summaries])...)
    if keep_original
        return CausalTable(merge(Tables.columns(x), sumtbl), x.treatment, x.response, x.controls, getgraph(x), getsummaries(x))
    else
        return sumtbl
    end
end

# Warning: this function is very slow and should be used only when there are no better options.
"""
    apply_function_over_neighbors(x::CausalTable, var_to_summarize::Symbol, summary_func::Function; use_inneighbors = true)

Apply a summary function over the neighbors of each row in a CausalTable. Mostly used as a utility to define a `summarize` function for new NetworkSummary types.

# Arguments
- `x::CausalTable`: The CausalTable object.
- `var_to_summarize::Symbol`: The symbol representing the variable in the CausalTable to summarize.
- `summary_func::Function`: The summary function to apply. Should take in a vector as input.
- `use_inneighbors::Bool`: (optional) If `true`, use the in-neighbors of each row as neighbors. If `false`, use the out-neighbors. Default is `true`.

# Returns
A Vector containing the summary value for each row.

"""
function apply_function_over_neighbors(x::CausalTable, var_to_summarize::Symbol, summary_func::Function; use_inneighbors = true)
    variable = Tables.getcolumn(x, var_to_summarize)
    if use_inneighbors
        return [summary_func(variable[inneighbors(x.graph, i)]) for i in 1:DataAPI.nrow(x)]
    else
        return [summary_func(variable[outneighbors(x.graph, i)]) for i in 1:DataAPI.nrow(x)]
    end
end

summarize(x::CausalTable, summary::NeighborSum) = adjacency_matrix(x.graph) * Tables.getcolumn(x, summary.var_to_summarize)
summarize(x::CausalTable, summary::NeighborProduct) = all(Tables.getcolumn(x, summary.var_to_summarize) .> 0) ? exp.(adjacency_matrix(x.graph) * log.(Tables.getcolumn(x, summary.var_to_summarize))) : apply_function_over_neighbors(x, summary.var_to_summarize, prod; use_inneighbors = summary.use_inneighbors)
summarize(x::CausalTable, summary::NeighborMaximum) = apply_function_over_neighbors(x, summary.var_to_summarize, maximum; use_inneighbors = summary.use_inneighbors)
summarize(x::CausalTable, summary::NeighborMinimum) = apply_function_over_neighbors(x, summary.var_to_summarize, minimum; use_inneighbors = summary.use_inneighbors)
summarize(x::CausalTable, summary::NeighborMode) = apply_function_over_neighbors(x, summary.var_to_summarize, StatsBase.mode; use_inneighbors = summary.use_inneighbors)
summarize(x::CausalTable, summary::Friends) = adjacency_matrix(x.graph) * ones(DataAPI.nrow(x))




