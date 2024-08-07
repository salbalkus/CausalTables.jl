"""
    abstract type NetworkSummary

Abstract type representing a summary of a network.

"""
abstract type NetworkSummary end
abstract type NetworkSummaryUnivariate <: NetworkSummary end
abstract type NetworkSummaryMultivariate <: NetworkSummary end

"""
    Sum <: NetworkSummary

A NetworkSummary which sums the values of the target variable for each unit connected in the adjacency `matrix`.

# Fields
- `target::Symbol`: The target variable of the network.
- `matrix::Symbol`: The matrix representation of the network.

"""
mutable struct Sum <: NetworkSummaryUnivariate
    target::Symbol
    matrix::Symbol
    weights::Union{Symbol, Nothing}
end
Sum(target::Symbol, matrix::Symbol) = Sum(target, matrix, nothing)
summarize(o::NamedTuple, x::Sum) = o[x.matrix] * (isnothing(x.weights) ? o[x.target] : o[x.target] .* o[x.weights])


mutable struct Mean <: NetworkSummaryUnivariate
    target::Symbol
    matrix::Symbol
    weights::Union{Symbol, Nothing}
end

Mean(target::Symbol, matrix::Symbol) = Mean(target, matrix, nothing)
function summarize(o::NamedTuple, x::Mean)
    denom = Base.replace(1 ./ vec(sum(o[x.matrix], dims = 2)), Inf => 0)
    v = o[x.matrix] * (isnothing(x.weights) ? o[x.target] : o[x.target] .* o[x.weights])
    return v .* denom
end

mutable struct AllOrderStatistics <: NetworkSummaryMultivariate
    target::Symbol
    matrix::Symbol
end

summarize(o::NamedTuple, x::AllOrderStatistics) = order_statistic_matrix(o[x.target], o[x.matrix])

function order_statistic_matrix(X::AbstractArray, G::SparseMatrixCSC)
    n = length(X)
    max_k = Int(maximum(sum(G, dims = 2)))
    Xvec = map(x -> Missings.missings(Float64, max_k), 1:n)
    r = rowvals(G)

    for i in 1:n
        ind = r[nzrange(G, i)]
        Xvec[i][1:length(ind)] = sort(X[ind], rev = true)
    end

    return(Xvec)
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
    summarize(o::CausalTable)

Summarizes the data in a `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object to be summarized.

# Returns
- A new `CausalTable` object with the original data merged with the summarized data.

"""
function summarize(o::CausalTable)
    
    # If we summarize a CausalTable that has already been summarized,
    # we need to replace the previous summary results.
    scm_result = merge(o.data, o.arrays)
    new_treatment = o.treatment
    new_confounders = o.confounders
    new_response = o.response

    nsummaries = length(o.summaries)
    tables = Vector{NamedTuple}(undef, nsummaries)

    for i in 1:nsummaries
        # get name and object of each summary
        sm = o.summaries[i]
        nm = keys(o.summaries)[i]

        # get the vector or matrix output by summary
        output_array = summarize(scm_result, sm)

        # construct a table, including header of names based on the summary name
        if eltype(output_array) <: AbstractArray
            header = map(x -> Symbol(string(nm, x)), (1:length(output_array[1])))
            tables[i] = Tables.columntable(Tables.table(mapreduce(permutedims, vcat, output_array); header = header))
        elseif ndims(output_array) == 1
            header = [nm]
            tables[i] = NamedTuple{(header...,)}((output_array,))
        else
            throw(ArgumentError("Summaries of arrays of dimension >2 not currently supported."))
        end

        # update the appropriate causal variable list based on the target 
        # of the summary with the newly generated columns

        origin = gettarget(sm)
        if origin ∈ o.treatment
            new_treatment = union(new_treatment, header)
        elseif origin ∈ o.confounders
            new_confounders = union(new_confounders, header)
        elseif origin ∈ o.response
            new_response = union(new_response, header)
        end
    end

    new_table = merge(o.data, tables...)

    return CausalTable(new_table, new_treatment, new_response, new_confounders, o.arrays, o.summaries)
end

gettarget(s::Friends) = nothing
gettarget(s::NetworkSummary) = s.target