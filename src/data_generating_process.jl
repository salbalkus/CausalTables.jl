
# Overload the convolve function to work on a vector of UnivariateDistribution
function Distributions.convolve(ds::Vector{T}) where {T <: UnivariateDistribution}
    output = ds[1]
    for d in ds[2:end]
        output = Distributions.convolve(output, d)
    end
    return output
end

### Summary functions ###
abstract type NetworkSummary end
abstract type NeighborSum <: NetworkSummary end

struct NeighborSumOut <: NeighborSum
    var_to_summarize::Symbol
end

struct NeighborSumIn <: NeighborSum
    var_to_summarize::Symbol
end

ValidDGPTypes = Union{Function, NetworkSummary}

### Constructors ###
mutable struct DataGeneratingProcess
    networkgen::Function
    distgen::Vector{Pair{Symbol, ValidDGPTypes}}
end

DataGeneratingProcess(distgen::Vector{Pair{Symbol, ValidDGPTypes}}) = DataGeneratingProcess(() -> nothing, distgen)

### Funtions to step through the DGP ###
# Step for summarizing a previous step over graph neighbors
initialize_dgp_step(step_func::NeighborSum, ct) = ct -> adjacency_matrix(ct.graph) * Tables.getcolumn(ct, step_func.var_to_summarize)

# Step for creating a distribution
initialize_dgp_step(step_func::Function, ct) = step_func(; ct.tbl...)

# Handle errors
function initialize_dgp_step(step_func, ct)
    error("Attempted to step through the DGP but failed due to incorrect type. Note that `step_func` must be of type Function or an existing NetworkSummary (e.g. NeighborSum)")
end

### Functions to draw a new column in the DGP ###
function append_dgp_draw!(name::Symbol, step::UnivariateDistribution, ct::CausalTable, n)
    ct.tbl = merge(ct.tbl, NamedTuple{(name,)}((Distributions.rand(step, n),)))
end

function append_dgp_draw!(name::Symbol, step::Vector{T}, ct::CausalTable, n) where {T <: UnivariateDistribution}
    ct.tbl = merge(ct.tbl, NamedTuple{(name,)}((map(Distributions.rand, step),)))
end

function append_dgp_draw!(name::Symbol, step::MultivariateDistribution, ct::CausalTable, n)
    ct.tbl = merge(ct.tbl, NamedTuple{(name,)}((Distributions.rand(step),)))
end

function append_dgp_draw!(name::Symbol, step::Function, ct::CausalTable, n) # append the summary function and store it in the CausalTable
    ct.tbl = merge(ct.tbl, NamedTuple{(name,)}((step(ct),)))
    ct.summaries = merge(ct.summaries, NamedTuple{(name,)}((step,)))
end

### Functions to get conditional distribution ###
get_conditional_distribution(varfunc::Function, ct) = varfunc(; ct.tbl...)

function get_conditional_distribution(varfunc::NeighborSumOut, ct)
    neighbordists = dgp.distgen[findfirst(isequal(varfunc.var_to_summarize), [name for (name, _) in dgp.distgen])][2](; ct.tbl...)
    return [Distributions.convolve(neighbordists[outneighbors(ct.graph, i)]) for i in 1:nv(ct.graph)] 
end

function get_conditional_distribution(varfunc::NeighborSumIn, ct)
    neighbordists = dgp.distgen[findfirst(isequal(varfunc.var_to_summarize), [name for (name, _) in dgp.distgen])][2](; ct.tbl...)
    return [Distributions.convolve(neighbordists[inneighbors(ct.graph, i)]) for i in 1:nv(ct.graph)] 
end

### Random Draw from DGP ###

function Base.rand(dgp::DataGeneratingProcess, n::Int)
    # Initialize output
    ct = CausalTable((;), dgp.networkgen(n), (;))

    # Iterate through each step of the DGP
    for pair in dgp.distgen
        print(pair)
        # Create the distribution (if step is a function)
        # or summarize a previous step over graph neighbors (if step is a NetworkSummary)
        dgp_step = initialize_dgp_step(pair[2], ct)

        # Draw from the distribution or pass the summarized data on
        append_dgp_draw!(pair[1], dgp_step, ct, n)
    end

    return ct
end

### Conditional Density of DGP ###

function condensity(dgp::DataGeneratingProcess, ct::CausalTable, var::Symbol)
    # only take the first instance of the variable name in the DGP
    varfunc = dgp.distgen[findfirst(isequal(var), [name for (name, _) in dgp.distgen])][2]
    return get_conditional_distribution(varfunc, ct)
end












