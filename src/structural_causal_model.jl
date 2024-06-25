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
    treatment::Vector{Symbol}
    response::Vector{Symbol}
    confounders::Vector{Symbol}

    function StructuralCausalModel(dgp, treatment, response, confounders)
        treatment, response, confounders = _process_causal_variable_names(treatment, response, confounders)
        new(dgp, treatment, response, confounders)
    end
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
    result = Vector{Pair{Symbol, Any}}(undef, length(scm))
    tag = Vector{Symbol}(undef, length(scm))
    summaries = (;)

    # Iterate through each step of the SCM
    for i_step in 1:length(scm)

        # Draw from the result of the step function
        output_step = scm.funcs[i_step](output...)
        result[i_step] = scm.dgp.names[i_step] => _scm_draw(output_step, output, n)
        
        # Determine where the resulting output should be placed in the CausalTable
        if result[i_step] <: AbstractVector{<:Number}
            tag[i_step] = :data
        elseif result[i_step] <: AbstractArray{<:Number, 2}
            tag[i_step] = :arrays
        else
            tag[i_step] = :nothing
        end

        # If the initial output was a NetworkSummary, record this to include in the CausalTable
        if output_step <: NetworkSummary
            summaries = merge(summaries, NamedTuple{scm.name[i_step]}((output_step,)))
        end    
    end

    data = NamedTuple(result[tag .== :data])
    arrays = NamedTuple(result[tag .== :arrays])

    # Store the recorded draws in a CausalTable format
    return CausalTable(data, scm.treatment, scm.response, scm.confounders, arrays, summaries)
end

# Helper functions for drawing a random CausalTable
_scm_draw(x::Distribution,            o::CausalTable, n::Int64) = rand(x, n)
_scm_draw(x::Vector{Distribution},    o::CausalTable, n::Int64) = length(x) == n ? rand.(x) : throw(ArgumentError("Length of vector of distributions in DataGeneratingProcess must be equal to n"))
_scm_draw(x::NetworkSummary,          o::CausalTable, n::Int64) = summarize(o, x)
_scm_draw(x::AbstractArray{<:Number}, o::CausalTable, n::Int64) = x


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
    varpos = findfirst(scm.names .== var)
    scm_result = getscm(ct)

    if scm.types[varpos] == :distribution
        return scm.funcs[varpos](scm_result)
    elseif scm.types[varpos] == :transformation
        return get_conditional_distribution(scm.funcs[varpos](scm_result), scm_result)
    else
        throw(ArgumentError("Cannot get conditional density. Variable $var is not a distribution or transformation of distributions in the StructuralCausalModel."))
    end
end

# Get the conditional density of a Sum in the DGP
function get_conditional_distribution(ns::Sum, scm::StructuralCausalModel, scm_result::NamedTuple)
    
    # Find the position of the target variable in the SCM
    targetpos = findfirst(scm.names .== ns.target)

    # Get the distribution of the target variable
    targetdist = scm.funcs[targetpos](scm_result)

    # Get the matrix involved in the sum
    m = scm_result[ns.matrix]

    # Compute the conditional distribution of the sum using convolution formula defined in utilities.jl
    return [Distributions.convolve(targetdist[m[row, :]]) for row in 1:size(m, 1)]
end

# Fallback for when no closed-form distribution is implemented
get_conditional_distribution(ns::T, scm_result::NamedTuple) where {T <: NetworkSummary} = throw(ArgumentError("No closed-form distribution is currently implemented for this NetworkSummary."))


"""
    conmean(dgp::DataGeneratingProcess, ct::CausalTable, var::Symbol)

Compute the conditional mean of a variable in a CausalTable based on a DataGeneratingProcess.

# Arguments
- `dgp::DataGeneratingProcess`: The DataGeneratingProcess object representing the data generating process.
- `ct::CausalTable`: The CausalTable object representing the data.
- `var::Symbol`: The variable for which to compute the conditional mean.

# Returns
An array of conditional means for the specified variable.

"""
conmean(dgp::DataGeneratingProcess, ct::CausalTable, var::Symbol) = mean.(condensity(dgp, ct, var))