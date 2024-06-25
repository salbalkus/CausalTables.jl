"""
    abstract type NetworkSummary

Abstract type representing a summary of a network.

"""
abstract type NetworkSummary end

mutable struct Sum <: NetworkSummary
    target::Symbol
    matrix::Symbol
end

summarize(o, x::Sum) = o.arrays[x.matrix] * o.data[x.target]