"""
    condensity(scm::StructuralCausalModel, ct::CausalTable, name::Symbol)

Compute the conditional density of variable `name` in CausalTable `ct` that has been drawn from StructuralCausalModel `scm`.

# Arguments
- `scm::StructuralCausalModel`: The StructuralCausalModel representing the data generating process.
- `ct::CausalTable`: The CausalTable containing the observed data.
- `name::Symbol`: The variable for which to compute the conditional density.

# Returns
The conditional density of the variable `var` given the observed data.

"""
function condensity(scm::StructuralCausalModel, ct::CausalTable, name::Symbol)
    
    varpos = findfirst(scm.dgp.names .== name)
    isnothing(varpos) && throw(ArgumentError("Variable $(name) is not contained within the StructuralCausalModel"))

    prev_names = scm.dgp.names[1:(varpos-1)]
    scm_result = NamedTupleTools.select(getscm(ct), prev_names)

    try
        if scm.dgp.types[varpos] == :distribution
            return scm.dgp.funcs[varpos](scm_result)
        elseif scm.dgp.types[varpos] == :transformation
            return get_conditional_distribution(scm.dgp.funcs[varpos](scm_result), scm, scm_result)
        else
            throw(ArgumentError("Cannot get conditional density. Variable $(name) is not a distribution or transformation of distributions in the StructuralCausalModel."))
        end
    catch e
        error(DIST_ERR_MSG(name))
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
    return convolve(targetdist, m)
end

# Fallback for when no closed-form distribution is implemented
get_conditional_distribution(ns::T, scm::StructuralCausalModel, scm_result::NamedTuple) where {T <: NetworkSummary} = throw(ArgumentError("No closed-form distribution is currently implemented for this NetworkSummary."))

"""
    conmean(scm::StructuralCausalModel, ct::CausalTable, name::Symbol)

Compute the conditional mean of variable `name` in CausalTable `ct` that has been drawn from StructuralCausalModel `scm`.

# Arguments
- `scm::StructuralCausalModel`: The StructuralCausalModel object representing the data generating process.
- `ct::CausalTable`: The CausalTable object representing the data.
- `name::Symbol`: The variable for which to compute the conditional mean.

# Returns
An array of conditional means for the specified variable.

"""
conmean(scm::StructuralCausalModel, ct::CausalTable, name::Symbol) = mean.(condensity(scm, ct, name))

"""
    convar(scm::StructuralCausalModel, ct::CausalTable, name::Symbol)

Compute the conditional variance of variable `name` in CausalTable `ct` that has been drawn from StructuralCausalModel `scm`.

# Arguments
- `scm::StructuralCausalModel`: The StructuralCausalModel object representing the data generating process.
- `ct::CausalTable`: The CausalTable object representing the data.
- `name::Symbol`: The variable for which to compute the conditional mean.

# Returns
An array of conditional variances for the specified variable.

"""
convar(scm::StructuralCausalModel, ct::CausalTable, name::Symbol) = var.(condensity(scm, ct, name))

"""
    propensity(scm::StructuralCausalModel, ct::CausalTable, name::Symbol)

Compute the (generalized) propensity score of variable `name` in CausalTable `ct` that has been drawn from StructuralCausalModel `scm`.

# Arguments
- `scm::StructuralCausalModel`: The StructuralCausalModel object representing the data generating process.
- `ct::CausalTable`: The CausalTable object representing the data.
- `name::Symbol`: The variable for which to compute the conditional mean.

# Returns
An array of conditional probabilities for the specified variable (or densities, if the specified variable is continuous).

"""
propensity(scm::StructuralCausalModel, ct::CausalTable, name::Symbol) = pdf.(condensity(scm, ct, name), Tables.getcolumn(ct, name))




