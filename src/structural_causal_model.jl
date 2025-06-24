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

    # Construct a Vector as subsequent inputs into each step function of the DGP
    result = (;)
    tag = Vector{Symbol}(undef, length(scm))
    summaries = (;)

    # Iterate through each step of the SCM
    for i_step in 1:length(scm)
        # Draw from the result of the step function
        result_step = try
            scm.dgp.funcs[i_step](result)
        catch e
            error(DIST_ERR_MSG(scm.dgp.names[i_step]))
        end
        result_draw =  CausalTables._scm_draw(result_step, result, n)
        result = merge(result, NamedTuple{(scm.dgp.names[i_step],)}((result_draw,)))
            
        # Determine where the resulting output should be placed in the CausalTable
        typeof_result_step = typeof(result_step)
        if (typeof_result_step <: Distribution) || (typeof_result_step <: AbstractArray{<:Distribution})
            tag[i_step] = :data
        else
            tag[i_step] = :arrays
        end
        # If the initial output was a NetworkSummary, record this to include in the CausalTable
        if typeof(result_step) <: CausalTables.NetworkSummary
            summaries = merge(summaries, NamedTuple{(scm.dgp.names[i_step],)}((result_step,)))
        end    
    end

    # Collect the variables to be stored as covariates in the CausalTable
    # We need to split apart any that generate matrices (e.g. draws from MvNormal)
    # from those that generate vectors (e.g. draws from Normal)
    data_tag = (tag .== :data)
    data_values = values(result)[data_tag]
    data_names = keys(result)[data_tag]

    # Iterate through the data, either re-storing the vectors or converting matrices 
    # to column tables, which are NamedTuples with each column as a separate variable,
    # and merging them into a single NamedTuple to store in the CausalTable
    data_output = (;)
    for (name, val) in zip(data_names, data_values)
        if length(size(val)) == 1
            # If the variable is a vector, we can just store it directly
            data_output = merge(data_output, NamedTuple{(name,)}((val,)))
        else
            new_val = Tables.table(val; header = [Symbol(name, i) for i in 1:size(val, 2)])
            # Store each column as a separate variable in the output
            data_output = NamedTupleTools.merge(data_output, Tables.columntable(new_val))
        end
    end

    # Collect extra arrays
    array_tag = (tag .!= :data)
    arrays = NamedTuple{keys(result)[array_tag]}(values(result)[array_tag],)    

    # Store the recorded draws in a CausalTable format
    return CausalTable(data_output, scm.treatment, scm.response, scm.causes, arrays, summaries)
end

Base.rand(scm::StructuralCausalModel) = rand(scm, 1)

### Helper functions for drawing a random CausalTable ###
# Single Distributions: Draw n samples
_scm_draw(x::T, o::NamedTuple, n::Int64) where {T <: UnivariateDistribution}     = rand(x, n)
_scm_draw(x::T, o::NamedTuple, n::Int64) where {T <: MultivariateDistribution}   = transpose(rand(x, n))
_scm_draw(x::T, o::NamedTuple, n::Int64) where {T <: MatrixDistribution}         = rand(x)

# Vectors of Distributions: Draw a single sample (assumes user input sample n into the distribution)
function _scm_draw(x::AbstractArray{<:Distribution}, o::NamedTuple, n::Int64)
        return rand.(x)
end

# Summary Function: Compute the summary
_scm_draw(x::T, o::NamedTuple, n::Int64) where {T<:NetworkSummary} = summarize(o, x)

# Fallback: Pass previous result directly
_scm_draw(x, o::NamedTuple, n::Int64) = x

"""
    update_arrays(scm::StructuralCausalModel, ct::CausalTable)

Propagate updates (e.g. interventions on treatment) through the `:code` portions of a Structural Causal Model (SCM), generating a new `CausalTable`.

# Arguments
- `scm::StructuralCausalModel`: The Structural Causal Model describing how each variable was drawn.
- `ct::CausalTable`: The `CausalTable` to update

# Returns
A `CausalTable` object containing the generated data.

"""
function update_arrays(scm::StructuralCausalModel, ct::CausalTable)

    scm_result = (;)
    new_arrays = (;)

    # Iterate through each step of the SCM
    for i_step in 1:length(scm)
        type = scm.dgp.types[i_step]
        name = scm.dgp.names[i_step]
        if type == :code
            scm_update = scm.dgp.funcs[i_step](scm_result)
            new_arrays = merge(new_arrays, NamedTuple{(name,)}((scm_update,)))
        elseif type == :distribution
            scm_update = Tables.getcolumn(ct, name)
        elseif type == :transformation
            scm_update = summarize(scm_result, scm.dgp.funcs[i_step](scm_result))
            new_arrays = merge(new_arrays, NamedTuple{(name,)}((scm_update,)))
        end
        scm_result = merge(scm_result, NamedTuple{(name,)}((scm_update,)))
    end

    # Store the recorded draws in a CausalTable format
    return CausalTable(ct.data, ct.treatment, ct.response, ct.causes, new_arrays, ct.summaries)
end



