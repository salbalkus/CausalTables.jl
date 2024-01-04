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

"""
struct NeighborSumIn <: NeighborSum
    var_to_summarize::Symbol
end

summarize(x::CausalTable) = (;zip(propertynames(x.summaries),  [summarize(x, s) for s in x.summaries])...)
summarize(x::CausalTable, summary::NeighborSum) = adjacency_matrix(x.graph) * Tables.getcolumn(x, summary.var_to_summarize)
