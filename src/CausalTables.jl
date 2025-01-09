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
include("conditional_density.jl")
include("estimands.jl")


# Exports

# causal_table.jl
export CausalTable, istable, columns, columnaccess
export nrow, ncol, subset, select, reject
export replace, treatment, response, confounders, data
export treatmentmatrix, responsematrix, treatmentparents, responseparents, parents
export confounders, confoundernames, confoundersmatrix
export mediators, mediatornames, mediatorsmatrix
export instruments, instrumentnames, instrumentsmatrix
export adjacency_matrix, dependency_matrix

# network_summary.jl
export NetworkSummary, Sum, Friends, Mean, AllOrderStatistics, KOrderStatistics
export summarize

# data_generating_process.jl
export DataGeneratingProcess
export @dgp

# structural_causal_model.jl
export StructuralCausalModel, getscm
export rand, condensity, conmean, convar, propensity

# estimands.jl
export intervene, draw_counterfactual 
export additive_mtp, multiplicative_mtp
export cfmean, cfdiff, ate, att, atu, ape
export treat_all, treat_none

end
