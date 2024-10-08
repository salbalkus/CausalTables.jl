module CausalTables

# Imports
using Distributions
using Tables
using TableTransforms
using PrettyTables
using DataAPI
using StatsBase
using LinearAlgebra
using SparseArrays
using Missings
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
include("estimands.jl")


# Exports

# causal_table.jl
export CausalTable, istable, columns, columnaccess
export nrow, ncol, subset
export replace, treatment, response, confounders
export treatmentnames, responsenames, confoundernames 
export treatmentparents, responseparents

# network_summary.jl
export NetworkSummary, Sum, Friends, Mean, AllOrderStatistics, KOrderStatistics
export summarize

# data_generating_process.jl
export DataGeneratingProcess
export @dgp

# structural_causal_model.jl
export StructuralCausalModel, getscm
export rand, condensity, conmean, convar

# estimands.jl
export draw_counterfactual, additive_mtp, multiplicative_mtp
export cfmean, cfdiff, ate, att, atu, ape
export treat_all, treat_none, cast_matrix_to_table_function

end
