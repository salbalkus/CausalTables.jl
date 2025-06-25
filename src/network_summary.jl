"""
    abstract type NetworkSummary

Abstract type representing a summary of a network.

"""
abstract type NetworkSummary end
abstract type NetworkSummaryUnivariate <: NetworkSummary end
abstract type NetworkSummaryMultivariate <: NetworkSummary end

"""
    Sum <: NetworkSummary

A NetworkSummary which sums the values of the target variable for each unit connected in the adjacency `matrix` of a `StructuralCausalModel` or `CausalTable`

# Fields
- `target::Symbol`: A key denoting the target variable to be summarized in the `DataGeneratingProcess` of the `StructuralCausalModel`; or, alternatively, the target variable to be summarized in the `data` attribute of a `CausalTable`.
- `matrix::Symbol`: A key denoting the adjacency matrix over which summary is computed in the `DataGeneratingProcess` of the `StructuralCausalModel`; or, alternatively, the key of the adjacency matrix in the `arrays` attribute of a `CausalTable`.
- `weights::Union{Symbol, Nothing}`: An optional variable by which each unit may be weighted in the summary.

"""
mutable struct Sum <: NetworkSummaryUnivariate
    target::Symbol
    matrix::Symbol
    weights::Union{Symbol, Nothing}
end

Sum(target::Symbol, matrix::Symbol) = Sum(target, matrix, nothing)
summarize(o::NamedTuple, x::Sum) = o[x.matrix] * (isnothing(x.weights) ? o[x.target] : o[x.target] .* o[x.weights])

"""
    Mean <: NetworkSummary

A NetworkSummary which computes the mean of the target variable among each unit connected in the adjacency `matrix`.

# Fields
- `target::Symbol`: A key denoting the target variable to be summarized in the `DataGeneratingProcess` of the `StructuralCausalModel`; or, alternatively, the target variable to be summarized in the `data` attribute of a `CausalTable`.
- `matrix::Symbol`: A key denoting the adjacency matrix over which summary is computed in the `DataGeneratingProcess` of the `StructuralCausalModel`; or, alternatively, the key of the adjacency matrix in the `arrays` attribute of a `CausalTable`.
- `weights::Union{Symbol, Nothing}`: An optional variable by which each unit may be weighted in the summary.

"""
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

"""
    AllOrderStatistics <: NetworkSummary

A NetworkSummary which computes all ordered values of the target variable among each unit's connected neighbors in the adjacency `matrix`.

# Fields
- `target::Symbol`: A key denoting the target variable to be summarized in the `DataGeneratingProcess` of the `StructuralCausalModel`; or, alternatively, the target variable to be summarized in the `data` attribute of a `CausalTable`.
- `matrix::Symbol`: A key denoting the adjacency matrix over which summary is computed in the `DataGeneratingProcess` of the `StructuralCausalModel`; or, alternatively, the key of the adjacency matrix in the `arrays` attribute of a `CausalTable`.
- `weights::Union{Symbol, Nothing}`: An optional variable by which each unit may be weighted in the summary.

"""
mutable struct AllOrderStatistics <: NetworkSummaryMultivariate
    target::Symbol
    matrix::Symbol
end

summarize(o::NamedTuple, x::AllOrderStatistics) = order_statistic_matrix(o[x.target], o[x.matrix] .!= 0, Int(maximum(sum(o[x.matrix] .!= 0, dims = 2))), true)

"""
    KOrderStatistics <: NetworkSummary

A NetworkSummary which computes the top K ordered values of the target variable among each unit's connected neighbors in the adjacency `matrix`.

# Fields
- `target::Symbol`: A key denoting the target variable to be summarized in the `DataGeneratingProcess` of the `StructuralCausalModel`; or, alternatively, the target variable to be summarized in the `data` attribute of a `CausalTable`.
- `matrix::Symbol`: A key denoting the adjacency matrix over which summary is computed in the `DataGeneratingProcess` of the `StructuralCausalModel`; or, alternatively, the key of the adjacency matrix in the `arrays` attribute of a `CausalTable`.
- `weights::Union{Symbol, Nothing}`: An optional variable by which each unit may be weighted in the summary.

"""
mutable struct KOrderStatistics <: NetworkSummaryMultivariate
    target::Symbol
    matrix::Symbol
    K::Int
    is_maximum::Bool
end
KOrderStatistics(target, matrix, K) = KOrderStatistics(target, matrix, K, true)

summarize(o::NamedTuple, x::KOrderStatistics) = order_statistic_matrix(o[x.target], o[x.matrix], x.K, x.is_maximum)

function order_statistic_matrix(X::AbstractArray, G::SparseMatrixCSC, max_k::Int, is_maximum::Bool)
    n = length(X)
    Xvec = map(x -> Missings.missings(Float64, max_k), 1:n)
    r = rowvals(G)

    for i in 1:n
        ind = r[nzrange(G, i)]
        cur_Xvec = sort(X[ind], rev = is_maximum)

        Xvec[i][1:minimum([length(cur_Xvec), max_k])] = (length(cur_Xvec) ≤ max_k ? cur_Xvec : cur_Xvec[1:max_k])
    end

    return(Xvec)
end

"""
    mutable struct Friends <: NetworkSummary

A NetworkSummary counting the number of connected individuals in an adjacency matrix, also known as the number of "friends".

# Fields
- `matrix::Symbol`: A key denoting the adjacency matrix over which summary is computed in the `DataGeneratingProcess` of the `StructuralCausalModel`; or, alternatively, the key of the adjacency matrix in the `arrays` attribute of a `CausalTable`.

"""
mutable struct Friends <: NetworkSummary
    matrix::Symbol
end

summarize(o::NamedTuple, x::Friends) = o[x.matrix] * ones(size(o[x.matrix], 1))



"""
    summarize(o::CausalTable)

Summarizes the data in a `CausalTable` object according to the NetworkSummary objects stored in its `summaries` attribute.

# Arguments
- `o::CausalTable`: The `CausalTable` object to be summarized.

# Returns
- A new `CausalTable` object with the original data merged with the summarized data.

"""
function summarize(o::CausalTable; add_summaries_as_causes = false)
    
    # If we summarize a CausalTable that has already been summarized,
    # we need to replace the previous summary results.
    scm_result = merge(o.data, o.arrays)
    new_treatment = o.treatment
    new_causes = deepcopy(o.causes)
    new_response = o.response

    nsummaries = length(o.summaries)
    tables = Vector{NamedTuple}(undef, nsummaries)
    headers = Vector{Vector}(undef, nsummaries)

    for i in 1:nsummaries
        # get name and object of each summary
        sm = o.summaries[i]
        nm = keys(o.summaries)[i]

        # get the vector or matrix output by summary
        output_array = summarize(scm_result, sm)

        # construct a table, including header of names based on the summary name
        if eltype(output_array) <: AbstractArray
            header = map(x -> Symbol(string(nm, x)), (1:length(output_array[1])))
            tables[i] = Tables.columntable(Tables.table(transpose(reduce(hcat, output_array)); header = header))
        elseif ndims(output_array) == 1
            header = [nm]
            tables[i] = NamedTuple{(header...,)}((output_array,))
        else
            throw(ArgumentError("Summaries of arrays of dimension >2 not currently supported."))
        end

        # update the appropriate causal variable list based on the target 
        # of the summary with the newly generated columns

        origin = CausalTables.gettarget(sm)
        if origin ∈ o.treatment
            new_treatment = union(new_treatment, header)
        elseif origin ∈ o.response
            new_response = union(new_response, header)
        end

        # update cause tuple
        if origin ∈ keys(o.causes)
            new_causes = merge(new_causes, NamedTuple{Tuple(header)}(repeat([deepcopy(o.causes[origin])], length(header))))
        end
        headers[i] = header
    end

    # setup to update causes by adding summarized symbols next to original symbols
    # need to extract headers from iteration above because some summaries add multiple new variables
    targets = CausalTables.gettarget.(values(o.summaries))
    originals = vcat([repeat([targets[i]], length(headers[i])) for i in 1:length(headers) if !isnothing(targets[i])]...)
    additions = vcat([headers[i] for i in 1:length(headers) if !isnothing(targets[i])]...)

    # option to add each summary as a cause if its original variables was also a cause
    if add_summaries_as_causes
        for cause in new_causes
            for c in cause
                append!(cause, additions[c .== originals])
            end
        end
    end

    # create a new data table and construct CausalTable
    new_table = merge(o.data, tables...)
    return CausalTable(new_table, new_treatment, new_response, new_causes, o.arrays, o.summaries)
end

gettarget(s::Friends) = nothing
gettarget(s::NetworkSummary) = s.target