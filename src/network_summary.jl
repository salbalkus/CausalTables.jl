"""
    abstract type NetworkSummary

Abstract type representing a summary of a network.

"""
abstract type NetworkSummary end

get_var_to_summarize(x::NetworkSummary) = x.var_to_summarize

"""
    abstract type NeighborSum <: NetworkSummary

Struct representing a summary of neighbor variables.

"""
mutable struct NeighborSum <: NetworkSummary 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    NeighborSum(var_to_summarize::Symbol; use_inneighbors::Bool = true) = new(var_to_summarize, use_inneighbors)
end

mutable struct NeighborProduct <: NetworkSummary 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    NeighborProduct(var_to_summarize::Symbol; use_inneighbors::Bool = true) = new(var_to_summarize, use_inneighbors)
end

abstract type NeighborOrderStatistic <: NetworkSummary end

mutable struct NeighborMaximum <: NeighborOrderStatistic 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    NeighborMaximum(var_to_summarize::Symbol; use_inneighbors::Bool = true) = new(var_to_summarize, use_inneighbors)
end

mutable struct NeighborMinimum <: NeighborOrderStatistic 
    var_to_summarize::Symbol
    use_inneighbors::Bool
    NeighborMinimum(var_to_summarize::Symbol; use_inneighbors::Bool = true) = new(var_to_summarize, use_inneighbors)
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

summarize(x::CausalTable, summary::NeighborSum) = adjacency_matrix(x.graph) * Tables.getcolumn(x, summary.var_to_summarize)

# TODO: Not sure if these work. Need to test
summarize(x::CausalTable, summary::NeighborProduct) = exp.(adjacency_matrix(x.graph) * log.(Tables.getcolumn(x, summary.var_to_summarize)))
summarize(x::CausalTable, summary::NeighborMaximum) = maximum(adjacency_matrix(x.graph) .* Tables.getcolumn(x, summary.var_to_summarize); dims = 2)
summarize(x::CausalTable, summary::NeighborMinimum) = minimum(adjacency_matrix(x.graph) .* Tables.getcolumn(x, summary.var_to_summarize); dims = 2)





