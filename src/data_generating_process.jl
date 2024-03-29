NET_ERR_MSG = "Failed to generate a valid network. This is likely because the function provided to the DataGeneratingProcess 
                does not return a valid Graph object. Ensure that `networkgen` is a function of the form n -> Graph(n, ...) 
                where `n` is the number of nodes in the graph."
DIST_ERR_MSG(node) = "Failed to generate a valid distribution for variable $(node). This is likely because the function provided to the 
                      DataGeneratingProcess does not return a valid distribution object. If using the @dgp macro, ensure that the 
                        right hand side of `$(node) ~ Distribution` returns either a valid distribution or a vector of distributions.
                        Otherwise, ensure that the distribution functions are of the form (; O...) -> Distribution(...) where `O` is a 
                        NamedTuple of the observed variables. Finally, ensure that when referencing previous variables in the DGP, 
                        that they are preceded by a colon, as in `:X`"
SUPPORTED_SUMMARIES = "Sum"
"""
    macro dgp(args...)

Defines a macro for creating a data generating process (DGP) in CausalTables.

# Examples
```jldoctest
julia> distseq = dgp(
    L1 ~ DiscreteUniform(1, 5),
    A ~ (@. Normal(:L1, 1)),
    Y ~ (@. Normal(:A + 0.2 * :L1, 1))
)
julia> dgp = DataGeneratingProcess(distseq);
```

"""
macro dgp(args...)
    return Vector{Pair{Symbol, CausalTables.ValidDGPTypes}}([_parse_tilde(arg) for arg in args])
end

# Helper function to parse each line in the dgp macro
function _parse_tilde(expr)
    # If expr is of the form `var1 = NetworkSummary(var2)`, turn it into a Pair directly
    if expr.head == :(=) && (typeof(eval(expr.args[2])) <: NetworkSummary)
        return expr.args[1] => eval(expr.args[2])
    # If expr is of the form `var ~ distribution`, turn it into an anonymous function that returns said distribution
    elseif expr.args[1] == :~
        pair = expr.args[2] => eval(:((; O...) -> $(_parse_distribution(expr.args[3]))))
        if !(pair[1] isa String || pair[1] isa Symbol)
           throw(ArgumentError("Invalid variable name. Variable name must be a string or symbol."))
        end
        return pair
    else
        throw(ArgumentError("Invalid expression. Must be of the form `var ~ distribution` or `var2 = NetworkSummary(var1)`."))
    end
end

# Helper function to parse `~ Distribution` type phrases in the dgp macro into anonymous functions
function _parse_distribution(expr)
    return postwalk(expr) do s
        typeof(s)==QuoteNode && return (:(O[$s]))
        s
    end
end


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
    treatment::SymbolOrNothing
    response::SymbolOrNothing
    controls::VectorOfSymbolsOrNothing
    function DataGeneratingProcess(networkgen::Function, distgen::Vector{Pair{Symbol, T}}, treatment::SymbolOrNothing, response::SymbolOrNothing, controls::VectorOfSymbolsOrNothing) where {T <: ValidDGPTypes}
        varnames = [name for (name, _) in distgen] 
        if !isnothing(controls) && (treatment ∈ controls || response ∈ controls)
            throw(ArgumentError("Treatment and/or response cannot be the same as controls."))
        elseif (!isnothing(treatment) && treatment ∉ varnames) || (!isnothing(response) && response ∉ varnames) || (!isnothing(controls) && any([c ∉ varnames for c in controls]))
            throw(ArgumentError("Treatment and/or response names not found in distribution generators."))
        end
        
        return new(networkgen, Vector{Pair{Symbol, ValidDGPTypes}}(distgen), treatment, response, controls)
    end
end

# Constructors
function DataGeneratingProcess(networkgen::Function, distgen::Vector{Pair{Symbol, T}}; treatment::SymbolOrNothing = nothing, response::SymbolOrNothing = nothing, controls::VectorOfSymbolsOrNothing = nothing) where {T <: ValidDGPTypes}
    return DataGeneratingProcess(networkgen, distgen, treatment, response, controls)
end

DataGeneratingProcess(distgen::Vector{Pair{Symbol, T}}; treatment::SymbolOrNothing = nothing, response::SymbolOrNothing = nothing, controls::VectorOfSymbolsOrNothing = nothing) where {T <: ValidDGPTypes} = DataGeneratingProcess(n -> nothing, distgen; treatment = treatment, response = response, controls = controls)

# Getters
gettreatmentsymbol(dgp::DataGeneratingProcess) = dgp.treatment
getresponsesymbol(dgp::DataGeneratingProcess) = dgp.response
getcontrolssymbol(dgp::DataGeneratingProcess) = dgp.controls

# Helper function to initialize the DGP step depending on its type using multiple dispatch
# Note that the step_func from a DataGeneratingProcess must be either a NetworkSummary or Function (otherwise it won't construct);
# this is an invariant, so we don't need to check for other types
_initialize_dgp_step(label, step_func::NetworkSummary, ct) = step_func
_initialize_dgp_step(label, step_func::Function, ct) = step_func(; ct.tbl...)


# Helper functions for drawing from the distribution or summarizing the data in a CausalTable using the `rand` function via multiple dispatch
function _append_dgp_draw!(name::Symbol, step::UnivariateDistribution, ct::CausalTable, n)
    ct.tbl = merge(ct.tbl, NamedTuple{(name,)}((Distributions.rand(step, n),)))
end
function _append_dgp_draw!(name::Symbol, step::Vector{T}, ct::CausalTable, n) where {T <: UnivariateDistribution}
    ct.tbl = merge(ct.tbl, NamedTuple{(name,)}((map(Distributions.rand, step),)))
end
function _append_dgp_draw!(name::Symbol, step::MultivariateDistribution, ct::CausalTable, n)
    ct.tbl = merge(ct.tbl, NamedTuple{(name,)}((Distributions.rand(step),)))
end
function _append_dgp_draw!(name::Symbol, step::NetworkSummary, ct::CausalTable, n) # summarize
    ct.tbl = merge(ct.tbl, NamedTuple{(name,)}((summarize(ct, step),)))
    ct.summaries = merge(ct.summaries, NamedTuple{(name,)}((step,)))
end

# Helper function to get the conditional distribution of a variable in a CausalTable depending on whether it is a distribution or a summary using multiple dispatch
_get_conditional_distribution(label, varfunc::Function, dgp::DataGeneratingProcess, ct::CausalTable) = try 
    varfunc(; ct.tbl...)
catch e
    error(DIST_ERR_MSG(label))
end

function _get_conditional_distribution(label, varfunc::Sum, dgp::DataGeneratingProcess, ct::CausalTable)
    try 
        neighbordists = dgp.distgen[findfirst(isequal(varfunc.var_to_summarize), [name for (name, _) in dgp.distgen])][2](; ct.tbl...)

        if varfunc.use_inneighbors
            return [
                indegree(ct.graph, i) > 0 ? Distributions.convolve(neighbordists[Graphs.inneighbors(ct.graph, i, varfunc.include_self)]) : Binomial(0) 
                for i in 1:nv(ct.graph) 
                ]
        else
            return [
                outdegree(ct.graph, i) > 0 ? Distributions.convolve(neighbordists[Graphs.outneighbors(ct.graph, i, varfunc.include_self)]) : Binomial(0) 
                for i in 1:nv(ct.graph) 
                ]
        end
    catch
        error("Sum cannot be computed over the specified distribution of $(label), likely because there is no closed-form convolution formula. Ensure that the distribution of $(label) adheres to the list given here: https://juliastats.org/Distributions.jl/stable/convolution/")
    end
end

# Fallbacks in case of error
function _get_conditional_distribution(label, varfunc::T, dgp::DataGeneratingProcess, ct::CausalTable) where {T <: NetworkSummary}
    error("NetworkSummary type not supported for computation over conditional distribution. Currently supported types are: $(SUPPORTED_SUMMARIES)")
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
    net = try
        dgp.networkgen(n)
    catch e
        error(NET_ERR_MSG)
    end
    !(typeof(net) <: Union{AbstractGraph, Nothing}) && error(NET_ERR_MSG)

    if (typeof(net) <: AbstractGraph)
        ct = CausalTable((;), gettreatmentsymbol(dgp), getresponsesymbol(dgp), getcontrolssymbol(dgp), net, (;))
    else
        ct = CausalTable((;), gettreatmentsymbol(dgp), getresponsesymbol(dgp), getcontrolssymbol(dgp), Graph(), (;))
    end

    # Iterate through each step of the DGP
    for pair in dgp.distgen
        # Create the distribution (if step is a function) or pass on the summary function
        dgp_step = try
            _initialize_dgp_step(pair[1], pair[2], ct)
        catch
            error(DIST_ERR_MSG(pair[1]))
        end

        # Draw from the distribution or summarize the data
        _append_dgp_draw!(pair[1], dgp_step, ct, n)
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
    label, varfunc = dgp.distgen[findfirst(isequal(var), [name for (name, _) in dgp.distgen])]
    return _get_conditional_distribution(label, varfunc, dgp, ct)
end

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













