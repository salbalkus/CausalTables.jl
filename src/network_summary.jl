"""
    abstract type NetworkSummary

Abstract type representing a summary of a network.

"""
abstract type NetworkSummary end

"""
    Sum <: NetworkSummary

A NetworkSummary which sums the values of the target variable for each unit connected in the adjacency `matrix`.

# Fields
- `target::Symbol`: The target variable of the network.
- `matrix::Symbol`: The matrix representation of the network.

"""
mutable struct Sum <: NetworkSummary
    target::Symbol
    matrix::Symbol
    weights::Union{Symbol, Nothing}
end
Sum(target::Symbol, matrix::Symbol) = Sum(target, matrix, nothing)
summarize(o::NamedTuple, x::Sum) = o[x.matrix] * (isnothing(x.weights) ? o[x.target] : o[x.target] .* o[x.weights])


mutable struct Mean <: NetworkSummary
    target::Symbol
    matrix::Symbol
    weights::Union{Symbol, Nothing}
end

Mean(target::Symbol, matrix::Symbol) = Mean(target, matrix, nothing)
function summarize(o::NamedTuple, x::Mean)
    denom = Base.replace(vec(sum(o[x.matrix], dims = 2)), Inf => 0)
    v = o[x.matrix] * (isnothing(x.weights) ? o[x.target] : o[x.target] .* o[x.weights])
    return v ./ denom
end


"""
    mutable struct Friends <: NetworkSummary

A NetworkSummary counting the number of connected individuals in an adjacency matrix, also known as the number of "friends"

# Fields
- `matrix::Symbol`: The matrix representing the network.

"""
mutable struct Friends <: NetworkSummary
    matrix::Symbol
end

summarize(o::NamedTuple, x::Friends) = o[x.matrix] * ones(size(o[x.matrix], 1))



"""
    summarize(o::CausalTable; keep_original=true)

Summarizes the data in a `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object to be summarized.
- `keep_original::Bool`: Whether to keep the original data in the resulting `CausalTable`. Default is `true`.

# Returns
- If `keep_original` is `true`, a new `CausalTable` object with the original data merged with the summarized data.
- If `keep_original` is `false`, a dictionary containing the summarized data.

"""
function summarize(o::CausalTable; keep_original = true)
    scm_result = merge(o.data, o.arrays)
    sum_data = (;zip(propertynames(o.summaries),  [summarize(scm_result, s) for s in o.summaries])...)
    if keep_original
        return replace(o, data = merge(o.data, sum_data))
    else
        return sum_data
    end
end

gettarget(s::Friends) = nothing
gettarget(s::NetworkSummary) = s.target