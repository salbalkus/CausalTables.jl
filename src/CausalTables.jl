module CausalTables

# Imports
using Distributions
using Graphs
using Tables
import Base: getindex

# Includes
include("utilities.jl")
include("causal_table.jl")
include("network_summary.jl")

ValidDGPTypes = Union{Function, NetworkSummary}

include("data_generating_process.jl")



# Exports

# causal_table.jl
export CausalTable, istable, columns, columnaccess, getcolumn, columnnames, columnindex, columntype
export getindex

# data_generating_process.jl
export NetworkSummary, NeighborSum, NeighborSumOut, NeighborSumIn
export DataGeneratingProcess, rand, condensity

end
