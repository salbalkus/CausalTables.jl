"""
    abstract type NetworkSummary

Abstract type representing a summary of a network.

"""
abstract type NetworkSummary end

"""
    abstract type NeighborSum <: NetworkSummary

Abstract type representing a summary of neighbor variables.

"""
abstract type NeighborSum <: NetworkSummary end

"""
    struct NeighborSumOut <: NeighborSum

Struct representing a summary of outgoing neighbor variables.

# Fields
- `var_to_summarize::Symbol`: The variable to summarize.

"""
struct NeighborSumOut <: NeighborSum
    var_to_summarize::Symbol
end

"""
    struct NeighborSumIn <: NeighborSum

Struct representing a summary of incoming neighbor variables.

# Fields
- `var_to_summarize::Symbol`: The variable to summarize.
- `keep::Bool`: Whether to append the summarized columns to the existing CausalTable. Defaults to `true`.

"""
struct NeighborSumIn <: NeighborSum
    var_to_summarize::Symbol
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
        return CausalTable(merge(Tables.columns(x), sumtbl), x.treatment, x.response, getgraph(x), getsummaries(x))
    else
        return sumtbl
    end
end

summarize(x::CausalTable, summary::NeighborSum) = adjacency_matrix(x.graph) * Tables.getcolumn(x, summary.var_to_summarize)
