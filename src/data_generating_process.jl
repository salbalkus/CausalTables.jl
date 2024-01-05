
"""
    struct DataGeneratingProcess

A mutable struct representing a data generating process.

# Fields
- `networkgen::Function`: A function that generates the network structure.
- `distgen::Vector{Pair{Symbol, ValidDGPTypes}}`: A vector of pairs representing the distribution generators for each variable.
- `treatment::Symbol`: The symbol representing the treatment variable.
- `response::Symbol`: The symbol representing the response variable.
- `controls::Vector{Symbol}`: A vector of symbols representing the control variables.

"""
mutable struct DataGeneratingProcess
    networkgen::Function
    distgen::Vector{Pair{Symbol, ValidDGPTypes}}
    treatment::Union{Symbol, Nothing}
    response::Union{Symbol, Nothing}
    controls::Union{Vector{Symbol}, Nothing}
end

# Constructors
DataGeneratingProcess(distgen::Vector{Pair{Symbol, T}}, treatment::Symbol, response::Symbol, controls::Vector{Symbol}) where {T <: ValidDGPTypes} = DataGeneratingProcess(n -> nothing, distgen, treatment, response, controls)
DataGeneratingProcess(distgen::Vector{Pair{Symbol, T}}) where {T <: ValidDGPTypes} = DataGeneratingProcess(n -> nothing, distgen, nothing, nothing, nothing)
DataGeneratingProcess(networkgen::Function, distgen::Vector{Pair{Symbol, T}}) where {T <: ValidDGPTypes} = DataGeneratingProcess(networkgen, distgen, nothing, nothing, nothing)

function DataGeneratingProcess(networkgen::Function, distgen::Vector{Pair{Symbol, T}}; treatment = nothing, response = nothing, controls = nothing) where {T <: ValidDGPTypes}
    if !isnothing(controls) && (treatment ∈ controls || response ∈ controls)
        error("Treatment and/or response cannot be the same as controls.")
    end
    return DataGeneratingProcess(networkgen, distgen, treatment, response, controls)
end

# Getters
gettreatment(dgp::DataGeneratingProcess) = dgp.treatment
getresponse(dgp::DataGeneratingProcess) = dgp.response
getcontrols(dgp::DataGeneratingProcess) = dgp.controls

"""
    initialize_dgp_step(step_func::NeighborSum, ct)

Initialize the data generating process (DGP) step by applying the given `step_func` to the causal table `ct`.

# Arguments
- `step_func`: The function to be applied to the causal table. Either a function or a NetworkSummary (e.g. NeighborSum).
- `ct`: The causal table.

# Returns
- `ct`: The modified causal table after applying the `step_func`.

"""
initialize_dgp_step(step_func::NetworkSummary, ct) = step_func
initialize_dgp_step(step_func::Function, ct) = step_func(; ct.tbl...)

function initialize_dgp_step(step_func, ct)
    error("Attempted to step through the DGP but failed due to incorrect type. Note that `step_func` must be of type Function or an existing NetworkSummary (e.g. NeighborSum)")
end

"""
    append_dgp_draw!(name::Symbol, step::UnivariateDistribution, ct::CausalTable, n)

Appends a draw from a distribution to the CausalTable.

# Arguments
- `name::Symbol`: The name of the column to be added.
- `step: Either the distribution (or a Vector of distributions) from which to draw values, or a function that summarizes neighboring distributions in the causal graph between units.
- `ct::CausalTable`: The CausalTable to which the values will be appended.
- `n`: The number of values to draw.

"""
function append_dgp_draw!(name::Symbol, step::UnivariateDistribution, ct::CausalTable, n)
    ct.tbl = merge(ct.tbl, NamedTuple{(name,)}((Distributions.rand(step, n),)))
end

function append_dgp_draw!(name::Symbol, step::Vector{T}, ct::CausalTable, n) where {T <: UnivariateDistribution}
    ct.tbl = merge(ct.tbl, NamedTuple{(name,)}((map(Distributions.rand, step),)))
end

function append_dgp_draw!(name::Symbol, step::MultivariateDistribution, ct::CausalTable, n)
    ct.tbl = merge(ct.tbl, NamedTuple{(name,)}((Distributions.rand(step),)))
end

function append_dgp_draw!(name::Symbol, step::NetworkSummary, ct::CausalTable, n) # summarize
    ct.tbl = merge(ct.tbl, NamedTuple{(name,)}((summarize(ct, step),)))
    ct.summaries = merge(ct.summaries, NamedTuple{(name,)}((step,)))
end

"""
    get_conditional_distribution(varfunc::Function, ct)

Compute the conditional distribution of a variable in the DataGeneratingProcess using a given function and a CausalTable.

# Arguments
- `varfunc`: The right side of the Pair used to compute the conditional distribution in the DataGeneratingProcess. Either of type Function or NetworkSummary (e.g. NeighborSum).
- `ct::CausalTable`: The CausalTable containing the data.

# Returns
The conditional distribution of the variable.

"""
get_conditional_distribution(varfunc::Function, dgp::DataGeneratingProcess, ct::CausalTable) = varfunc(; ct.tbl...)

function get_conditional_distribution(varfunc::NeighborSumOut, dgp::DataGeneratingProcess, ct::CausalTable)
    neighbordists = dgp.distgen[findfirst(isequal(varfunc.var_to_summarize), [name for (name, _) in dgp.distgen])][2](; ct.tbl...)
    return [
        outdegree(ct.graph, i) > 0 ? Distributions.convolve(neighbordists[outneighbors(ct.graph, i)]) : Binomial(0) 
        for i in 1:nv(ct.graph) 
        ] end

function get_conditional_distribution(varfunc::NeighborSumIn, dgp::DataGeneratingProcess, ct::CausalTable)
    neighbordists = dgp.distgen[findfirst(isequal(varfunc.var_to_summarize), [name for (name, _) in dgp.distgen])][2](; ct.tbl...)
    return [
        indegree(ct.graph, i) > 0 ? Distributions.convolve(neighbordists[inneighbors(ct.graph, i)]) : Binomial(0) 
        for i in 1:nv(ct.graph) 
        ] 
end

"""
    rand(dgp::DataGeneratingProcess, n::Int)

Generate a random CausalTable using the specified DataGeneratingProcess.

# Arguments
- `dgp::DataGeneratingProcess`: The DataGeneratingProcess object defining the causal network and distribution steps.
- `n::Int`: The number of observations to generate.

# Returns
- `ct::CausalTable`: The generated CausalTable.

"""
function Base.rand(dgp::DataGeneratingProcess, n::Int)
    # Initialize output
    net = dgp.networkgen(n)
    if net isa SimpleGraph
        ct = CausalTable((;), gettreatment(dgp), getresponse(dgp), getcontrols(dgp), net, (;))
    else
        ct = CausalTable((;), gettreatment(dgp), getresponse(dgp), getcontrols(dgp), Graph(), (;))
    end

    # Iterate through each step of the DGP
    for pair in dgp.distgen
        # Create the distribution (if step is a function) or pass on the summary function
        dgp_step = initialize_dgp_step(pair[2], ct)

        # Draw from the distribution or summarize the data
        append_dgp_draw!(pair[1], dgp_step, ct, n)
    end

    return ct
end


"""
    condensity(dgp::DataGeneratingProcess, ct::CausalTable, var::Symbol)

Compute the conditional density of a variable in a CausalTable based on a DataGeneratingProcess.

# Arguments
- `dgp::DataGeneratingProcess`: The DataGeneratingProcess object representing the data generating process.
- `ct::CausalTable`: The CausalTable object containing the data.
- `var::Symbol`: The name of the variable for which to compute the conditional density.

# Returns
The conditional density of the variable in the CausalTable.

"""
function condensity(dgp::DataGeneratingProcess, ct::CausalTable, var::Symbol)
    # only take the first instance of the variable name in the DGP
    varfunc = dgp.distgen[findfirst(isequal(var), [name for (name, _) in dgp.distgen])][2]
    return get_conditional_distribution(varfunc, dgp, ct)
end












