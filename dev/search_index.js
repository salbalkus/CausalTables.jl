var documenterSearchIndex = {"docs":
[{"location":"#CausalTables.jl","page":"CausalTables.jl","title":"CausalTables.jl","text":"","category":"section"},{"location":"","page":"CausalTables.jl","title":"CausalTables.jl","text":"The following is the documentation for each function in the CausalTables.jl package.","category":"page"},{"location":"","page":"CausalTables.jl","title":"CausalTables.jl","text":"Modules = [CausalTables]\nOrder   = [:function, :type]","category":"page"},{"location":"#Base.rand-Tuple{DataGeneratingProcess, Int64}","page":"CausalTables.jl","title":"Base.rand","text":"rand(dgp::DataGeneratingProcess, n::Int)\n\nGenerate a random CausalTable using the specified DataGeneratingProcess.\n\nArguments\n\ndgp::DataGeneratingProcess: The DataGeneratingProcess object defining the causal network and distribution steps.\nn::Int: The number of observations to generate.\n\nReturns\n\nct::CausalTable: The generated CausalTable.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.append_dgp_draw!-Tuple{Symbol, Distributions.UnivariateDistribution, CausalTable, Any}","page":"CausalTables.jl","title":"CausalTables.append_dgp_draw!","text":"append_dgp_draw!(name::Symbol, step::UnivariateDistribution, ct::CausalTable, n)\n\nAppends a draw from a distribution to the CausalTable.\n\nArguments\n\nname::Symbol: The name of the column to be added.\n`step: Either the distribution (or a Vector of distributions) from which to draw values, or a function that summarizes neighboring distributions in the causal graph between units.\nct::CausalTable: The CausalTable to which the values will be appended.\nn: The number of values to draw.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.condensity-Tuple{DataGeneratingProcess, CausalTable, Symbol}","page":"CausalTables.jl","title":"CausalTables.condensity","text":"condensity(dgp::DataGeneratingProcess, ct::CausalTable, var::Symbol)\n\nCompute the conditional density of a variable in a CausalTable based on a DataGeneratingProcess.\n\nArguments\n\ndgp::DataGeneratingProcess: The DataGeneratingProcess object representing the data generating process.\nct::CausalTable: The CausalTable object containing the data.\nvar::Symbol: The name of the variable for which to compute the conditional density.\n\nReturns\n\nThe conditional density of the variable in the CausalTable.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.get_conditional_distribution-Tuple{Function, DataGeneratingProcess, CausalTable}","page":"CausalTables.jl","title":"CausalTables.get_conditional_distribution","text":"get_conditional_distribution(varfunc::Function, ct)\n\nCompute the conditional distribution of a variable in the DataGeneratingProcess using a given function and a CausalTable.\n\nArguments\n\nvarfunc: The right side of the Pair used to compute the conditional distribution in the DataGeneratingProcess. Either of type Function or NetworkSummary (e.g. NeighborSum).\nct::CausalTable: The CausalTable containing the data.\n\nReturns\n\nThe conditional distribution of the variable.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.getgraph-Tuple{CausalTable}","page":"CausalTables.jl","title":"CausalTables.getgraph","text":"getgraph(x::CausalTable)\n\nGet the graph associated with a CausalTable.\n\nArguments\n\nx::CausalTable: The CausalTable object.\n\nReturns\n\nThe graph associated with the CausalTable.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.getresponse-Tuple{CausalTable}","page":"CausalTables.jl","title":"CausalTables.getresponse","text":"getresponse(x::CausalTable)\n\nGet the response variable from a CausalTable.\n\nArguments\n\nx::CausalTable: The CausalTable object.\n\nReturns\n\nThe response variable from the CausalTable.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.getsummaries-Tuple{CausalTable}","page":"CausalTables.jl","title":"CausalTables.getsummaries","text":"getsummaries(x::CausalTable)\n\nReturns the tables stored in the CausalTable x.\n\nArguments\n\nx::CausalTable: The CausalTable object.\n\nReturns\n\nAn array of tables stored in the CausalTable x.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.gettreatment-Tuple{CausalTable}","page":"CausalTables.jl","title":"CausalTables.gettreatment","text":"gettreatment(x::CausalTable)\n\nGet the treatment column from a CausalTable.\n\nArguments\n\nx::CausalTable: The CausalTable object.\n\nReturns\n\nThe treatment column of the CausalTable.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.initialize_dgp_step-Tuple{NetworkSummary, Any}","page":"CausalTables.jl","title":"CausalTables.initialize_dgp_step","text":"initialize_dgp_step(step_func::NeighborSum, ct)\n\nInitialize the data generating process (DGP) step by applying the given step_func to the causal table ct.\n\nArguments\n\nstep_func: The function to be applied to the causal table. Either a function or a NetworkSummary (e.g. NeighborSum).\nct: The causal table.\n\nReturns\n\nct: The modified causal table after applying the step_func.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.summarize","page":"CausalTables.jl","title":"CausalTables.summarize","text":"summarize(x::CausalTable, keep = true)\n\nSummarize a CausalTable by merging its columns, treatment, response, graph, and summaries.\n\nArguments:\n\nx::CausalTable: The CausalTable to be summarized.\nkeep::Bool: Determines whether to keep the original CausalTable or return a new summarized CausalTable. Default is true.\n\nReturns:\n\nIf keep is true, a new CausalTable with merged columns, treatment, response, graph, and summaries.\nIf keep is false, a NamedTuple with summaries as keys and the corresponding summarized CausalTable as values.\n\n\n\n\n\n","category":"function"},{"location":"#Distributions.convolve-Union{Tuple{Vector{T}}, Tuple{T}} where T<:(Distributions.UnivariateDistribution)","page":"CausalTables.jl","title":"Distributions.convolve","text":"Distributions.convolve(ds::Vector{T}) where {T <: UnivariateDistribution}\n\nOverload the convolve function to work on a vector of UnivariateDistribution.\n\nArguments\n\nds::Vector{T}: A vector of UnivariateDistribution objects.\n\nReturns\n\noutput: The result of convolving all the distributions in ds.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.CausalTable","page":"CausalTables.jl","title":"CausalTables.CausalTable","text":"CausalTable\n\nA mutable structure that contains a table (tbl), a graph (graph), and a named tuple of summaries (summaries).\n\n\n\n\n\n","category":"type"},{"location":"#CausalTables.CausalTable-Tuple{Any, Graphs.SimpleGraphs.SimpleGraph, NamedTuple}","page":"CausalTables.jl","title":"CausalTables.CausalTable","text":"CausalTable(tbl)\n\nAn alternate constructor for CausalTable that takes a table, graph, and summaries, and initializes the treatment and response as blank.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.CausalTable-Tuple{Any, Union{Nothing, Symbol}, Union{Nothing, Symbol}}","page":"CausalTables.jl","title":"CausalTables.CausalTable","text":"CausalTable(tbl)\n\nAn alternate constructor for CausalTable that takes a table, response, and treatment, and initializes the graph and summaries as blank.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.CausalTable-Tuple{Any}","page":"CausalTables.jl","title":"CausalTables.CausalTable","text":"CausalTable(tbl)\n\nAn alternate constructor for CausalTable that takes a table (tbl) and initializes all other variables as blank.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.DataGeneratingProcess","page":"CausalTables.jl","title":"CausalTables.DataGeneratingProcess","text":"mutable struct DataGeneratingProcess\n\nMutable struct representing a data generating process.\n\nFields\n\nnetworkgen::Function: The function to generate the network. Should return a SimpleGraph object and take a single argument n::Int64, the number of nodes in the graph\ndistgen::Vector{Pair{Symbol, ValidDGPTypes}}: The vector of variable-summary pairs.\n\n\n\n\n\n","category":"type"},{"location":"#CausalTables.DataGeneratingProcess-Union{Tuple{Array{Pair{Symbol, T}, 1}}, Tuple{T}} where T<:Union{NetworkSummary, Function}","page":"CausalTables.jl","title":"CausalTables.DataGeneratingProcess","text":"DataGeneratingProcess(distgen::Vector{Pair{Symbol, ValidDGPTypes}})\n\nConstructs a DataGeneratingProcess object with the given distribution generator.\n\nArguments\n\ndistgen::Vector{Pair{Symbol, ValidDGPTypes}}: A vector of pairs representing the distribution generator.\n\nReturns\n\nDataGeneratingProcess: The constructed DataGeneratingProcess object.\n\n\n\n\n\n","category":"method"},{"location":"#CausalTables.NeighborSum","page":"CausalTables.jl","title":"CausalTables.NeighborSum","text":"abstract type NeighborSum <: NetworkSummary\n\nAbstract type representing a summary of neighbor variables.\n\n\n\n\n\n","category":"type"},{"location":"#CausalTables.NeighborSumIn","page":"CausalTables.jl","title":"CausalTables.NeighborSumIn","text":"struct NeighborSumIn <: NeighborSum\n\nStruct representing a summary of incoming neighbor variables.\n\nFields\n\nvar_to_summarize::Symbol: The variable to summarize.\nkeep::Bool: Whether to append the summarized columns to the existing CausalTable. Defaults to true.\n\n\n\n\n\n","category":"type"},{"location":"#CausalTables.NeighborSumOut","page":"CausalTables.jl","title":"CausalTables.NeighborSumOut","text":"struct NeighborSumOut <: NeighborSum\n\nStruct representing a summary of outgoing neighbor variables.\n\nFields\n\nvar_to_summarize::Symbol: The variable to summarize.\n\n\n\n\n\n","category":"type"},{"location":"#CausalTables.NetworkSummary","page":"CausalTables.jl","title":"CausalTables.NetworkSummary","text":"abstract type NetworkSummary\n\nAbstract type representing a summary of a network.\n\n\n\n\n\n","category":"type"}]
}
