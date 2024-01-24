"""
    abstract type NetworkSummary

Abstract type representing a summary of a network.

"""
abstract type NetworkSummary end

get_var_to_summarize(x::NetworkSummary) = x.var_to_summarize


"""
    mutable struct NeighborSum <: NetworkSummary

A mutable struct representing a network summary that calculates the sum of a variable's neighbors.

# Fields
- `var_to_summarize::Symbol`: The variable to summarize.
- `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors for the summary. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

# Constructors
- `NeighborSum(var_to_summarize::Symbol; use_inneighbors::Bool = true)`: Create a new NeighborSum object.

"""

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

"""
    mutable struct NeighborProduct <: NetworkSummary

A mutable struct representing a network summary that calculates the product of a variable's neighbors.

# Fields
- `var_to_summarize::Symbol`: The variable to summarize.
- `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

# Constructors
- `NeighborProduct(var_to_summarize::Symbol; use_inneighbors::Bool = true)`: Create a new `NeighborProduct` object.

"""
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

"""
    mutable struct NeighborMaximum <: NeighborOrderStatistic

A mutable struct representing a network summary that calculates the maximum of a variable's neighbors.

# Fields
- `var_to_summarize::Symbol`: The variable to summarize.
- `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

# Constructors
- `NeighborMaximum(var_to_summarize::Symbol; use_inneighbors::Bool = true)`: Create a new `NeighborMaximum` object.

"""
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

"""
    mutable struct NeighborMinimum <: NeighborOrderStatistic

A mutable struct representing a network summary that calculates the maximum of a variable's neighbors.

# Fields
- `var_to_summarize::Symbol`: The variable to summarize.
- `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

# Constructors
- `NeighborMinimum(var_to_summarize::Symbol; use_inneighbors::Bool = true)`: Create a new `NeighborMinimum` object.

"""
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

"""
    mutable struct NeighborMode <: NetworkSummary

A mutable struct representing a network summary that calculates the mode of a variable's neighbors.

# Fields
- `var_to_summarize::Symbol`: The variable to summarize.
- `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

# Constructors
- `NeighborMode(var_to_summarize::Symbol; use_inneighbors::Bool = true)`: Create a new `NeighborMode` object.

"""
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

"""
    mutable struct NeighborMode <: NetworkSummary

A mutable struct representing a network summary that calculates the number of neighbors of each unit.

# Fields
- `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

# Constructors
- `Friends(var_to_summarize::Symbol; use_inneighbors::Bool = true)`: Create a new `Friends` object.

"""
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

Summarize a CausalTable by merging its columns, treatment, response, graph, and summaries.

Arguments:
- `x::CausalTable`: The CausalTable to be summarized.
- `keep::Bool`: Determines whether to keep the original CausalTable or return a new summarized CausalTable. Default is `true`.

Returns:
- If `keep` is `true`, a new CausalTable with merged columns, treatment, response, graph, and summaries.
- If `keep` is `false`, a NamedTuple with summaries as keys and the corresponding summarized CausalTable as values.
"""
function summarize(x::CausalTable, keep = true)
    sumtbl = (;zip(propertynames(x.summaries),  [summarize(x, s) for s in x.summaries])...)
    if keep
        return CausalTable(merge(Tables.columns(x), sumtbl), x.treatment, x.response, x.controls, getgraph(x), getsummaries(x))
    else
        return sumtbl
    end
end

# Warning: this function is very slow and should be used only when there are no better options.
function apply_function_over_neighbors(x::CausalTable, var_to_summarize::Symbol, summary_func::Function; use_inneighbors = true)
    variable = Tables.getcolumn(x, var_to_summarize)
    if use_inneighbors
        return [summary_func(variable[inneighbors(x.graph, i)]) for i in 1:DataAPI.nrow(x)]
    else
        return [summary_func(variable[outneighbors(x.graph, i)]) for i in 1:DataAPI.nrow(x)]
    end
end

summarize(x::CausalTable, summary::NeighborSum) = adjacency_matrix(x.graph) * Tables.getcolumn(x, summary.var_to_summarize)
summarize(x::CausalTable, summary::NeighborProduct) = exp.(adjacency_matrix(x.graph) * log.(Tables.getcolumn(x, summary.var_to_summarize)))
summarize(x::CausalTable, summary::NeighborMaximum) = maximum(adjacency_matrix(x.graph) .* Tables.getcolumn(x, summary.var_to_summarize); dims = 2)
summarize(x::CausalTable, summary::NeighborMinimum) = minimum(adjacency_matrix(x.graph) .* Tables.getcolumn(x, summary.var_to_summarize); dims = 2)
summarize(x::CausalTable, summary::NeighborMode) = apply_function_over_neighbors(x, summary.var_to_summarize, StatsBase.mode; use_inneighbors = summary.use_inneighbors)
summarize(x::CausalTable, summary::Friends) = adjacency_matrix(x.graph) * ones(DataAPI.nrow(x))




