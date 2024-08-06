DIST_ERR_MSG(node) = "Failed to generate a valid distribution for variable $(node). This is likely because the function provided to the 
                      DataGeneratingProcess does not return a valid distribution object. If using the @dgp macro, ensure that the 
                        right hand side of `$(node) ~ Distribution` returns either a valid distribution or a vector of distributions.
                        Otherwise, ensure that the distribution functions are of the form O -> Distribution(...) where `O` is a 
                        NamedTuple of the observed variables."
SUPPORTED_SUMMARIES = "Sum"

# Helper function to check whether any elements of a vector are not contained in the DGP
not_in_dgp(dgp, vec) = any(map(x -> x ∉ dgp.names, vec))

"""
    struct StructuralCausalModel

A struct representing a structural causal model (SCM).

# Fields
- `dgp::DataGeneratingProcess`: The data generating process associated with the structural causal model.
- `treatment::Vector{Symbol}`: The variables representing the treatment in the structural causal model.
- `response::Vector{Symbol}`: The variables representing the response in the structural causal model.
- `confounders::Vector{Symbol}`: The variables representing the confounders in the structural causal model.

# Constructors
- `StructuralCausalModel(dgp, treatment, response, confounders)`: Constructs a new `StructuralCausalModel` object.

# Arguments
- `dgp`: The data generating process associated with the structural causal model.
- `treatment`: The variables representing the treatment in the structural causal model.
- `response`: The variables representing the response in the structural causal model.
- `confounders`: The variables representing the confounders in the structural causal model.
"""
mutable struct StructuralCausalModel
    dgp::DataGeneratingProcess
    treatment::Symbols
    response::Symbols
    confounders::Symbols
    arraynames::Symbols

    function StructuralCausalModel(dgp, treatment, response, confounders, arraynames)
        treatment, response, confounders = _process_causal_variable_names(treatment, response, confounders)

        not_in_dgp(dgp, treatment) && throw(ArgumentError("One or more of treatment labels $(treatment) are not a variable in the DataGeneratingProcess."))
        not_in_dgp(dgp, response) && throw(ArgumentError("One or more of response labels $(response) are not a variable in the DataGeneratingProcess."))
        not_in_dgp(dgp, confounders) && throw(ArgumentError("One or more of confounder labels $(confounders) are not a variable in the DataGeneratingProcess."))
        not_in_dgp(dgp, arraynames) && throw(ArgumentError("One or more of array labels $(arraynames) are not a variable in the DataGeneratingProcess."))

        new(dgp, treatment, response, confounders, arraynames)
    end
end

StructuralCausalModel(dgp, treatment, response, confounders; arraynames = []) = StructuralCausalModel(dgp, treatment, response, confounders, arraynames)
StructuralCausalModel(dgp, treatment, response; confounders = [], arraynames = []) = StructuralCausalModel(dgp, treatment, response, confounders, arraynames)

function StructuralCausalModel(dgp; treatment = nothing, response = nothing, confounders = [], arraynames = [])
    isnothing(treatment) && throw(ArgumentError("Treatment variable must be defined"))
    isnothing(response) && throw(ArgumentError("Response variable must be defined"))
    StructuralCausalModel(dgp, treatment, response, confounders, arraynames)
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

    # Collect the data
    data_tag = (tag .== :data)
    data = NamedTuple{keys(result)[data_tag]}(values(result)[data_tag],)

    # Collect extra arrays
    array_tag = (tag .!= :data)
    arrays = NamedTuple{keys(result)[array_tag]}(values(result)[array_tag],)    

    # Store the recorded draws in a CausalTable format
    return CausalTable(data, scm.treatment, scm.response, scm.confounders, arrays, summaries)
end

### Helper functions for drawing a random CausalTable ###
# Single Distributions: Draw n samples
_scm_draw(x::T, o::NamedTuple, n::Int64) where {T <: UnivariateDistribution}     = rand(x, n)
_scm_draw(x::T, o::NamedTuple, n::Int64) where {T <: MultivariateDistribution}   = rand(x)
_scm_draw(x::T, o::NamedTuple, n::Int64) where {T <: MatrixDistribution}         = rand(x)

# Vectors of Distributions: Draw a single sample (assumes user input sample n into the distribution)
function _scm_draw(x::AbstractArray{<:Distribution}, o::NamedTuple, n::Int64)
    if length(x) == n
        return rand.(x)
    else 
        throw(ArgumentError("Length of vector of distributions in DataGeneratingProcess must be equal to n"))
    end
end

# Summary Function: Compute the summary
_scm_draw(x::T, o::NamedTuple, n::Int64) where {T<:NetworkSummary} = summarize(o, x)

# Fallback: Pass previous result directly
_scm_draw(x, o::NamedTuple, n::Int64) = x


"""
    condensity(scm::StructuralCausalModel, ct::CausalTable, var::Symbol)

Compute the conditional density of a variable in a StructuralCausalModel given a CausalTable.

# Arguments
- `scm::StructuralCausalModel`: The StructuralCausalModel representing the data generating process.
- `ct::CausalTable`: The CausalTable containing the observed data.
- `var::Symbol`: The variable for which to compute the conditional density.

# Returns
The conditional density of the variable `var` given the observed data.

"""
function condensity(scm::StructuralCausalModel, ct::CausalTable, var::Symbol)
    
    varpos = findfirst(scm.dgp.names .== var)
    isnothing(varpos) && throw(ArgumentError("Argument `var` is not contained within the StructuralCausalModel"))
    
    scm_result = getscm(ct)

    try
        if scm.dgp.types[varpos] == :distribution
            return scm.dgp.funcs[varpos](scm_result)
        elseif scm.dgp.types[varpos] == :transformation
            return get_conditional_distribution(scm.dgp.funcs[varpos](scm_result), scm, scm_result)
        else
            throw(ArgumentError("Cannot get conditional density. Variable $var is not a distribution or transformation of distributions in the StructuralCausalModel."))
        end
    catch e
        error(DIST_ERR_MSG(var))
    end
end

# Get the conditional density of a Sum in the DGP
function get_conditional_distribution(ns::Sum, scm::StructuralCausalModel, scm_result::NamedTuple)
    
    # Find the position of the target variable in the SCM
    targetpos = findfirst(scm.dgp.names .== ns.target)

    # Get the distribution of the target variable
    targetdist = scm.dgp.funcs[targetpos](scm_result)

    # Cast the matrix involved in the sum to a Boolean matrix
    m = .!(iszero.(scm_result[ns.matrix]))

    # Compute the conditional distribution of the sum using convolution formula defined in utilities.jl
    return Distributions.convolve(targetdist, m)
end

# Fallback for when no closed-form distribution is implemented
get_conditional_distribution(ns::T, scm::StructuralCausalModel, scm_result::NamedTuple) where {T <: NetworkSummary} = throw(ArgumentError("No closed-form distribution is currently implemented for this NetworkSummary."))


"""
    conmean(scm::StructuralCausalModel, ct::CausalTable, var::Symbol)

Compute the conditional mean of a variable in a CausalTable based on a DataGeneratingProcess.

# Arguments
- `scm::StructuralCausalModel`: The StructuralCausalModel object representing the data generating process.
- `ct::CausalTable`: The CausalTable object representing the data.
- `var::Symbol`: The variable for which to compute the conditional mean.

# Returns
An array of conditional means for the specified variable.

"""
conmean(scm::StructuralCausalModel, ct::CausalTable, var::Symbol) = mean.(condensity(scm, ct, var))