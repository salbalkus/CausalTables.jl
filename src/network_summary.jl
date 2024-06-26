"""
    abstract type NetworkSummary

Abstract type representing a summary of a network.

"""
abstract type NetworkSummary end

get_var_to_summarize(x::NetworkSummary) = x.var_to_summarize

mutable struct Sum <: NetworkSummary 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    include_self::Bool
    
    @doc raw"""
        Sum(var_to_summarize::Symbol; use_inneighbors::Bool = true)

    Constructs a Network Summary denoting the **sum** of `var_to_summarize` over each unit's neighbors. If `var_to_summarize` is ``X``, then mathematically this computes

    ```math
    \sum_{j \in \mathcal{F}_i} X_j
    ```

    where ``F_i`` represents the "friends" or neighbors of unit ``i``.

    # Arguments
    - `var_to_summarize::Symbol`: The variable to summarize.
    - `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.
    - `include_self::Bool`: Whether to include the value of `var_to_summarize` for the unit itself in the sum. Default is `true`.

    # Returns
    - `Sum`: The constructed Sum object.
    """
    Sum(var_to_summarize::Symbol; use_inneighbors::Bool = true, include_self = true) = new(var_to_summarize, use_inneighbors, include_self)
end

mutable struct Product <: NetworkSummary 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    include_self::Bool

    @doc raw"""
        Product(var_to_summarize::Symbol; use_inneighbors::Bool = true)

    Constructs a Network Summary denoting the **product** of `var_to_summarize` over each unit's neighbors. If `var_to_summarize` is ``X``, then mathematically this computes

    ```math
    \prod_{j \in \mathcal{F}_i} X_j
    ```

    where ``F_i`` represents the "friends" or neighbors of unit ``i``.

    # Arguments
    - `var_to_summarize::Symbol`: The variable to summarize.
    - `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

    # Returns
    - `Product`: The constructed Product object.
    """
    Product(var_to_summarize::Symbol; use_inneighbors::Bool = true, include_self = true) = new(var_to_summarize, use_inneighbors, include_self)
end

abstract type OrderStatistic <: NetworkSummary end

mutable struct Maximum <: OrderStatistic 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    include_self::Bool

    @doc raw"""
        Maximum(var_to_summarize::Symbol; use_inneighbors::Bool = true)

    Constructs a Network Summary denoting the **maximum** of `var_to_summarize` over each unit's neighbors. If `var_to_summarize` is ``X``, then mathematically this computes

    ```math
    \max_{j \in \mathcal{F}_i} X_j
    ```

    where ``F_i`` represents the "friends" or neighbors of unit ``i``.

    # Arguments
    - `var_to_summarize::Symbol`: The variable to summarize.
    - `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

    # Returns
    - `Maximum`: The constructed Maximum object.
    """
    Maximum(var_to_summarize::Symbol; use_inneighbors::Bool = true, include_self::Bool = true) = new(var_to_summarize, use_inneighbors, include_self)
end

mutable struct Minimum <: OrderStatistic 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    include_self::Bool
    @doc raw"""
        Minimum(var_to_summarize::Symbol; use_inneighbors::Bool = true)

    Constructs a Network Summary denoting the **minimum** of `var_to_summarize` over each unit's neighbors. If `var_to_summarize` is ``X``, then mathematically this computes

    ```math
    \min_{j \in \mathcal{F}_i} X_j
    ```

    where ``F_i`` represents the "friends" or neighbors of unit ``i``.

    # Arguments
    - `var_to_summarize::Symbol`: The variable to summarize.
    - `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

    # Returns
    - `Minimum`: The constructed Minimum object.
    """
    Minimum(var_to_summarize::Symbol; use_inneighbors::Bool = true, include_self::Bool = true) = new(var_to_summarize, use_inneighbors, include_self)
end

mutable struct Mode <: NetworkSummary 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    include_self::Bool
    @doc raw"""
        Mode(var_to_summarize::Symbol; use_inneighbors::Bool = true)

    Constructs a Network Summary denoting the **mode** of `var_to_summarize` over each unit's neighbors -- that is, the value occuring most often among units connected to each unit in the network.

    # Arguments
    - `var_to_summarize::Symbol`: The variable to summarize.
    - `use_inneighbors::Bool`: Whether to use the in-neighbors or out-neighbors. If `true`, then the in-neighbors are used. If the graph is undirected, this has no effect.

    # Returns
    - `Mode`: The constructed Mode object.
    """
    Mode(var_to_summarize::Symbol; use_inneighbors::Bool = true, include_self::Bool = true) = new(var_to_summarize, use_inneighbors, include_self)
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
    summarize(x::CausalTable, keep_original = true)

Summarize a CausalTable by computing the summaries specified in the `summaries` field of the CausalTable and adding them to the existing data.

# Arguments
- `x::CausalTable`: The CausalTable to be summarized.
- `keep_original::Bool`: Whether to return a CausalTable with summarized variables appended (`true`) or a Tables-compatible NamedTuple with only the summarized variables (`false`).

# Returns
- If `keep_original` is `true`, a new CausalTable with merged columns, treatment, response, graph, and summaries. Otherwise, NamedTuple of just variable summaries

"""
function summarize(x::CausalTable; keep_original = true)
    sumtbl = (;zip(propertynames(x.summaries),  [summarize(x, s) for s in x.summaries])...)
    if keep_original
        return CausalTable(merge(Tables.columns(x), sumtbl), x.treatment, x.response, x.controls, getgraph(x), getsummaries(x))
    else
        return sumtbl
    end
end

# Overload neighbor functions from Graphs
Graphs.inneighbors(g::AbstractGraph, v::Int, include_self) = include_self ? [v; Graphs.inneighbors(g, v)] : Graphs.inneighbors(g, v)
Graphs.outneighbors(g::AbstractGraph, v::Int, include_self) = include_self ? [v; Graphs.outneighbors(g, v)] : Graphs.outneighbors(g, v)

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
function apply_function_over_neighbors(x::CausalTable, var_to_summarize::Symbol, summary_func::Function; use_inneighbors = true, include_self = true)
    variable = Tables.getcolumn(x, var_to_summarize)
    if use_inneighbors
        return [summary_func(variable[Graphs.inneighbors(x.graph, i, include_self)]) for i in 1:DataAPI.nrow(x)]
    else
        return [summary_func(variable[Graphs.outneighbors(x.graph, i, include_self)]) for i in 1:DataAPI.nrow(x)]
    end
end

summarize(x::CausalTable, summary::Sum) = try 
    adjacency_matrix(getgraph(x)) * Tables.getcolumn(x, summary.var_to_summarize) .+ (summary.include_self ? Tables.getcolumn(x, summary.var_to_summarize) : 0)
catch
    throw(ArgumentError("Attempted to summarize a variable over a graph that is not the correct size. In a CausalTable, the Graph must have the same number of vertices as the number of rows of the Table. 
    If no graph is being provided to the table, then it is not possible to apply a summary."))
end
function summarize(x::CausalTable, summary::Product)
    if all(Tables.getcolumn(x, summary.var_to_summarize) .> 0)
        output = exp.(adjacency_matrix(x.graph) * log.(Tables.getcolumn(x, summary.var_to_summarize)))
    else
        output = apply_function_over_neighbors(x, summary.var_to_summarize, prod; use_inneighbors = summary.use_inneighbors, include_self = summary.include_self)
    end
    return output
end
summarize(x::CausalTable, summary::Maximum) = apply_function_over_neighbors(x, summary.var_to_summarize, maximum; use_inneighbors = summary.use_inneighbors, include_self = summary.include_self)
summarize(x::CausalTable, summary::Minimum) = apply_function_over_neighbors(x, summary.var_to_summarize, minimum; use_inneighbors = summary.use_inneighbors, include_self = summary.include_self)
summarize(x::CausalTable, summary::Mode) = apply_function_over_neighbors(x, summary.var_to_summarize, StatsBase.mode; use_inneighbors = summary.use_inneighbors, include_self = summary.include_self)
summarize(x::CausalTable, summary::Friends) = adjacency_matrix(x.graph) * ones(DataAPI.nrow(x))


