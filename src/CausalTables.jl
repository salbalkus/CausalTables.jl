module CausalTables

# Imports
using Distributions
using Tables
using TableTransforms
using PrettyTables
using DataAPI
using StatsBase
using LinearAlgebra
import Base: getindex
import MacroTools: postwalk

# New types
Symbols = AbstractArray{Symbol, 1}

# Includes
include("utilities.jl")
include("causal_table.jl")
include("network_summary.jl")

include("data_generating_process.jl")
include("structural_causal_model.jl")

# Exports

# causal_table.jl
export CausalTable, istable, columns, columnaccess
export nrow, ncol, subset
export replace, treatment, response, confounders
export treatmentnames, responsenames, confoundernames 
export treatmentparents, responseparents

# network_summary.jl
export NetworkSummary, Sum, Friends
export summarize

# data_generating_process.jl
export DataGeneratingProcess
export @dgp

# structural_causal_model.jl
export StructuralCausalModel, getscm
export rand, condensity, conmean

end
