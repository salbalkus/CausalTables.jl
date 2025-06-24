DIST_ERR_MSG(node) = "Failed to generate a valid distribution for variable $(node). This is likely because the function provided to the DataGeneratingProcess does not return a valid distribution object. If using the @dgp macro, ensure that the right hand side of `$(node) ~ Distribution` returns either a valid distribution or a vector of distributions. Otherwise, ensure that the distribution functions are of the form O -> Distribution(...) where `O` is a NamedTuple of the observed variables."
SUPPORTED_SUMMARIES = "Sum"

# Helper function to check whether any elements of a vector are not contained in the DGP
not_in_dgp(dgp, vec) = any(map(x -> x ∉ dgp.names, vec))

"""
    struct StructuralCausalModel

A struct representing a structural causal model (SCM). This includes a DataGeneratingProcess 

# Arguments
- `dgp::DataGeneratingProcess`: The data generating process from which random data will be drawn.
- `treatment::Vector{Symbol}`: The variables representing the treatment.
- `response::Vector{Symbol}`: The variables representing the response.
- `causes::Union{NamedTuple, Nothing}`: A NamedTuple of Vectors labeling the causes of relevant variables in the data-generating process. If `nothing`, will assume that all variables not contained in `treatment` or `response` are common causes of both.
- `arraynames`: Names of auxiliary variables used in the DataGeneratingProcess that are not included as "tabular" variables. Most commonly used to denote names of adjacency matrices used to compute summary functions of previous steps. 

"""
mutable struct StructuralCausalModel
    dgp::DataGeneratingProcess
    treatment::Symbols
    response::Symbols
    causes::Union{NamedTuple, Nothing}
    arraynames::Symbols

    function StructuralCausalModel(dgp, treatment, response, causes, arraynames)
        treatment, response, causes = _process_causal_variable_names(treatment, response, causes)

        not_in_dgp(dgp, treatment) && throw(ArgumentError("One or more of treatment labels $(treatment) are not a variable in the DataGeneratingProcess."))
        not_in_dgp(dgp, response) && throw(ArgumentError("One or more of response labels $(response) are not a variable in the DataGeneratingProcess."))
        if(!isnothing(causes))
            not_in_dgp(dgp, union(keys(causes), vcat(values(causes)...))) && throw(ArgumentError("One or more symbols in `cause` are not a variable in the DataGeneratingProcess."))
        end
        not_in_dgp(dgp, arraynames) && throw(ArgumentError("One or more of array labels $(arraynames) are not a variable in the DataGeneratingProcess."))

        new(dgp, treatment, response, causes, arraynames)
    end
end

StructuralCausalModel(dgp, treatment::Symbol, response::Symbol, causes::NamedTuple; arraynames = []) = StructuralCausalModel(dgp, [treatment], [response], causes, arraynames)
StructuralCausalModel(dgp, treatment::Symbol, response::Symbols, causes::NamedTuple; arraynames = []) = StructuralCausalModel(dgp, [treatment], response, causes, arraynames)
StructuralCausalModel(dgp, treatment::Symbols, response::Symbol, causes::NamedTuple; arraynames = []) = StructuralCausalModel(dgp, treatment, [response], causes, arraynames)
StructuralCausalModel(dgp, treatment, response, causes; arraynames = []) = StructuralCausalModel(dgp, treatment, response, causes, arraynames)
StructuralCausalModel(dgp, treatment, response; causes = nothing, arraynames = []) = StructuralCausalModel(dgp, treatment, response, causes, arraynames)

function StructuralCausalModel(dgp; treatment = nothing, response = nothing, causes = nothing, arraynames = [])
    isnothing(treatment) && throw(ArgumentError("Treatment variable must be defined"))
    isnothing(response) && throw(ArgumentError("Response variable must be defined"))
    StructuralCausalModel(dgp, treatment, response, causes, arraynames)
end

Base.length(scm::StructuralCausalModel) = length(scm.dgp)

"""
    rand(scm::StructuralCausalModel, n::Int)

Generate random data from a Structural Causal Model (SCM) using the specified number of samples.

# Arguments
- `scm::StructuralCausalModel`: The Structural Causal Model from which to generate data.
- `n::Int`: The number of samples to generate.

# Returns
A `CausalTable` object containing the generated data.

"""
function Base.rand(scm::StructuralCausalModel, n::Int)

    ###  First, generate a direct sample from the DGP
    path = rand(scm.dgp, n)

    ### Second, process the path into a CausalTable
    ct_data = (;)

    # Iterate through each of the "objects" we've generated,
    # and decide where to place them in the CausalTable
    for i in 1:length(scm.dgp)
        cur_name = scm.dgp.names[i]

        # If the step in the path is a distribution,
        # put it in the data attribute of the CausalTable
        if scm.dgp.types[i] == :distribution
            if ndims(path[cur_name]) == 1
                ct_data = NamedTupleTools.merge(ct_data, NamedTuple{(cur_name,)}((path[cur_name],)))
            else
                # Split the matrix into a NamedTuple of vectors and merge it
                tbl = Tables.columntable(Tables.table(path[scm.dgp.names[i]]; header = [Symbol(scm.dgp.names[i], j) for j in 1:size(path[cur_name], 2)]))
                ct_data = NamedTupleTools.merge(ct_data, tbl)
            end
        end
    end

    # Put any parts of the path that were not already processed into the arrays attribute
    ct_arrays = NamedTupleTools.select(path, scm.dgp.names[scm.dgp.types .== :code])
    summaries = (;)

    return CausalTable(ct_data, scm.treatment, scm.response, scm.causes, ct_arrays, summaries)
end

Base.rand(scm::StructuralCausalModel) = rand(scm, 1)




