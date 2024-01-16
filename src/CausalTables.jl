module CausalTables

# Imports
using Distributions
using Graphs
using Tables
using TableOperations
using DataAPI
import Base: getindex
import MacroTools: postwalk

# Includes
include("utilities.jl")
include("causal_table.jl")
include("network_summary.jl")

ValidDGPTypes = Union{Function, NetworkSummary}
SymbolOrNothing = Union{Symbol, Nothing}
VectorOfSymbolsOrNothing = Union{Vector{Symbol}, Nothing}

include("data_generating_process.jl")

# Exports

# causal_table.jl
export CausalTable, istable, columns, columnaccess, getcolumn, columnnames, columnindex, columntype
export getindex, nrow, ncol
export getresponse, gettreatment, getsummaries, getgraph, getcontrols, gettable, getresponsesymbol, gettreatmentsymbol, getcontrolssymbols
export replacetable, replace, subset

# network_summary.jl
export NetworkSummary, NeighborSum
export summarize, get_var_to_summarize

# data_generating_process.jl
export DataGeneratingProcess, rand, condensity, conmean
export @scm

end
