"""
    abstract type NetworkSummary

Abstract type representing a summary of a network.

"""
abstract type NetworkSummary end

mutable struct Sum <: NetworkSummary
    target::Symbol
    matrix::Symbol
end

summarize(o::NamedTuple, x::Sum) = o[x.matrix] * o[x.target]


function summarize(o::CausalTable; keep_original = true)
    scm_result = merge(o.data, o.arrays)
    sum_data = (;zip(propertynames(o.summaries),  [summarize(scm_result, s) for s in o.summaries])...)
    if keep_original
        return replace(o, data = merge(o.data, sum_data))
    else
        return sum_data
    end
end