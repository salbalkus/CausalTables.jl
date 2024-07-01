var documenterSearchIndex = {"docs":
[{"location":"man/network-summaries/#Network-Summaries","page":"Network Summaries","title":"Network Summaries","text":"","category":"section"},{"location":"man/network-summaries/","page":"Network Summaries","title":"Network Summaries","text":"Typically, performing causal inference on a network relies on summarizing the treatment and covariates of each unit's neighbors using some sort of summary function. For example, in a study evaluating the effect of electric vehicle adoption on air pollution, one might model the commuting patterns between counties as a network and evaluate the effect of the sum the number of electric vehicles commuting into each county. The following documents summary measures available in CausalTables.jl and how to summarize a CausalTable","category":"page"},{"location":"man/network-summaries/#Summarizing-a-CausalTable","page":"Network Summaries","title":"Summarizing a CausalTable","text":"","category":"section"},{"location":"man/network-summaries/","page":"Network Summaries","title":"Network Summaries","text":"Data wrapped in a CausalTable includes a NamedTuple summaries which describes extra variables represented as summary variables over the network. These summary measures can be computed and added to the table by calling the summarize function on the CausalTable object.","category":"page"},{"location":"man/network-summaries/","page":"Network Summaries","title":"Network Summaries","text":"summarize","category":"page"},{"location":"man/network-summaries/#CausalTables.summarize","page":"Network Summaries","title":"CausalTables.summarize","text":"summarize(o::CausalTable; keep_original=true)\n\nSummarizes the data in a CausalTable object.\n\nArguments\n\no::CausalTable: The CausalTable object to be summarized.\nkeep_original::Bool: Whether to keep the original data in the resulting CausalTable. Default is true.\n\nReturns\n\nIf keep_original is true, a new CausalTable object with the original data merged with the summarized data.\nIf keep_original is false, a dictionary containing the summarized data.\n\n\n\n\n\n","category":"function"},{"location":"man/network-summaries/#Existing-Summary-Measures","page":"Network Summaries","title":"Existing Summary Measures","text":"","category":"section"},{"location":"man/network-summaries/","page":"Network Summaries","title":"Network Summaries","text":"The following lists summary measures currently available off-the-shelf in CausalTables.jl. Examples on their use are provided in Generating Data for Statistical Experiments and Turning data into a CausalTable.","category":"page"},{"location":"man/network-summaries/#Defining-Your-Own-Summary-Measures","page":"Network Summaries","title":"Defining Your Own Summary Measures","text":"","category":"section"},{"location":"man/network-summaries/","page":"Network Summaries","title":"Network Summaries","text":"Forthcoming.","category":"page"},{"location":"man/ground-truth/#Computing-Ground-Truth-Conditional-Distributions","page":"Computing Ground Truth of Causal Parameters","title":"Computing Ground Truth Conditional Distributions","text":"","category":"section"},{"location":"man/ground-truth/","page":"Computing Ground Truth of Causal Parameters","title":"Computing Ground Truth of Causal Parameters","text":"Once we've defined a DGP and have some table of data with variables matching those of our DGP, we can compute the ground truth conditional distributions of any variable in a CausalTable (given a corresponding DGP) using the condensity function. This returns a Distribution object from the package Distributions.jl.","category":"page"},{"location":"man/ground-truth/","page":"Computing Ground Truth of Causal Parameters","title":"Computing Ground Truth of Causal Parameters","text":"Let's see an example. First, we'll define a DGP:","category":"page"},{"location":"man/ground-truth/","page":"Computing Ground Truth of Causal Parameters","title":"Computing Ground Truth of Causal Parameters","text":"using Graphs\nusing CausalTables\nusing Random\nusing Distributions\n\ndgp = @dgp(\n        W ~ Binomial(10, 0.3),\n        X ~ (@. Normal(W + 1)),\n        A = adjacency_matrix(barabasi_albert(length(X), 2)),\n        Xs $ Sum(:X, :A),\n        Y ~ (@. LogNormal(log(0.2 * Xs + 4), 0.1 * W + 1))\n    )\n\nscm = StructuralCausalModel(\n    dgp;\n    treatment = :X,\n    response = :Y,\n    confounders = [:W]\n)","category":"page"},{"location":"man/ground-truth/","page":"Computing Ground Truth of Causal Parameters","title":"Computing Ground Truth of Causal Parameters","text":"Now, let's generate some data and compute the ground truth conditional distributions of the variables in the data. Note that if the DGP attempts to summarize a variable with no neighbors in a graph, the resulting conditional distribution will currently be Binomial(0, 0.5), which denotes a point-mass distribution at 0.","category":"page"},{"location":"man/ground-truth/","page":"Computing Ground Truth of Causal Parameters","title":"Computing Ground Truth of Causal Parameters","text":"Random.seed!(1);\ndata = rand(scm, 5)\nW_distribution = condensity(scm, data, :W)\nX_distribution = condensity(scm, data, :X)\nXs_distribution = condensity(scm, data, :Xs)","category":"page"},{"location":"man/ground-truth/","page":"Computing Ground Truth of Causal Parameters","title":"Computing Ground Truth of Causal Parameters","text":"One can also compute the ground truth conditional mean of a variable in a CausalTable using the conmean function:","category":"page"},{"location":"man/ground-truth/","page":"Computing Ground Truth of Causal Parameters","title":"Computing Ground Truth of Causal Parameters","text":"Y = conmean(scm, data, :Y)","category":"page"},{"location":"man/generating-data/#Generating-Data-for-Statistical-Experiments","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"","category":"section"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"When evaluating a causal inference method, we often want to test it on data from a known causal model. CausalTables.jl allows us to define a DataGeneratingProcess (or DGP) to do just that. ","category":"page"},{"location":"man/generating-data/#Defining-a-DataGeneratingProcess","page":"Generating Data for Statistical Experiments","title":"Defining a DataGeneratingProcess","text":"","category":"section"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"A data generating process describes a mechanism by which draws from random variables are simulated. It typically takes the form of a sequence of conditional distributions. CausalTables allows us to define a DGP as a DataGeneratingProcess object, which takes three arguments: the names of variables generated at each step, the types of these variables, and funcs, an array of functions of the form (; O...) -> *some code. ","category":"page"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"Suppose, for example, that we wanted to simulate data from the following DGP:","category":"page"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"beginalign*\n    W sim textDiscreteUniform(1 5) \n    X sim textNormal(W 1) \n    Y sim textNormal(X + 02W 1)\nendalign*","category":"page"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"where X is the treatment, Y is the response, and W is a confounding variable affecting both X and Y. A verbose and inconvenient (albeit correct) way to define this DGP would be as follows:","category":"page"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"using Distributions\nusing CausalTables\n\nDataGeneratingProcess(\n    [:W, :X, :Y],\n    [:distribution, :distribution, :distribution],\n    [\n        (; O...) -> DiscreteUniform(1, 5), \n        (; O...) -> (@. Normal(O.W, 1)),\n        (; O...) -> (@. Normal(O.X + 0.2 * O.W, 1))\n    ]\n)","category":"page"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"where ; O... syntax is a shorthand for a function that takes keyword arguments corresponding to the names of the variables in the DGP. ","category":"page"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"However, a much more convenient way to define this DGP is using the @dgp macro, which takes a sequence of conditional distributions of the form [variable name] ~ Distribution(args...) and deterministic variable assignments of the form [variable name] = f(...) and automatically generates a valid DataGeneratingProcess. For example, the easier way to define the DGP above is as follows:","category":"page"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"using CausalTables\ndistributions = @dgp(\n        W ~ DiscreteUniform(1, 5),\n        X ~ (@. Normal(W, 1)),\n        Y ~ (@. Normal(X + 0.2 * W, 1))\n    )","category":"page"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"Note that with the @dgp macro, any symbol (that is, any string of characters prefixed by a colon, as in :W or :X) is automatically replaced with the corresponding previously-defined variable in the process. For instance, in Normal(:W, 1), the :W will be replaced automatically with the distribution we defined as W earlier in the sequence. ","category":"page"},{"location":"man/generating-data/#Defining-a-StructuralCausalModel","page":"Generating Data for Statistical Experiments","title":"Defining a StructuralCausalModel","text":"","category":"section"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"In CausalTables.jl, a StructuralCausalModel is a data generating process endowed with some causal interpretation. Constructing a StructuralCausalModel allows users to randomly draw a CausalTable with the necessary components from the DataGeneratingProcess they've defined. With the above DataGeneratingProcess in hand, we can define a StructuralCausalModel object like so – treatment, response, and confounder variables in the causal model are specified as keyword arguments to the DataGeneratingProcess constructor:","category":"page"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"dgp = StructuralCausalModel(\n    distributions;\n    treatment = :X,\n    response = :Y,\n    confounders = [:W]\n)","category":"page"},{"location":"man/generating-data/#Networks-of-Causally-Connected-Units","page":"Generating Data for Statistical Experiments","title":"Networks of Causally-Connected Units","text":"","category":"section"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"In some cases, we might work with data in which units may not be causally independent, but rather, in which one unit's variables could dependent on some summary function of its neighbors. Generating data from such a model can be done by adding lines of the form Xs $ NetworkSummary to the @dgp macro.","category":"page"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"Here's an example of how such a DataGeneratingProcess might be constructed:","category":"page"},{"location":"man/generating-data/","page":"Generating Data for Statistical Experiments","title":"Generating Data for Statistical Experiments","text":"using Graphs\nusing CausalTables\nusing Distributions\n\ndgp = @dgp(\n        W ~ DiscreteUniform(1, 5),\n        n = length(W),\n        A = adjacency_matrix(erdos_renyi(n, 0.5)),\n        Ws $ Sum(:W, :A),\n        X ~ (@. Normal(Ws, 1)),\n        Xs $ Sum(:X, :A),\n        Y ~ (@. Normal(Xs + 0.2 * Ws, 1))\n    )\n\nscm = StructuralCausalModel(\n    dgp;\n    treatment = :X,\n    response = :Y,\n    confounders = [:W, :Ws]\n)","category":"page"},{"location":"man/formatting/#Turning-Your-Data-Into-a-CausalTable","page":"Turning data into a CausalTable","title":"Turning Your Data Into a CausalTable","text":"","category":"section"},{"location":"man/formatting/","page":"Turning data into a CausalTable","title":"Turning data into a CausalTable","text":"One of the main purposes of CausalTables.jl is to wrap a Table of data in Julia in order to provide it as input to some other causal inference package. Given a Table of some data, we can turn it into a CausalTable by specifying the treatment, response, and control variables. ","category":"page"},{"location":"man/formatting/#Tables-with-Causally-Independent-Units","page":"Turning data into a CausalTable","title":"Tables with Causally Independent Units","text":"","category":"section"},{"location":"man/formatting/","page":"Turning data into a CausalTable","title":"Turning data into a CausalTable","text":"The code below demonstrates this on the Titanic dataset. This could be, for example, to use as input into some estimator of whether a passenger's sex caused them to survive the Titanic disaster, controlling for some baselineline confounders listed in confounders.","category":"page"},{"location":"man/formatting/","page":"Turning data into a CausalTable","title":"Turning data into a CausalTable","text":"using CausalTables\nusing MLDatasets: Titanic\nusing DataFrames\n\ndf = Titanic().dataframe\n\n# Wrapping the dataset in a CausalTable\nctbl = CausalTable(df; treatment = :Sex, response = :Survived, confounders = [:Pclass, :Age, :SibSp])\n\nnothing # hide","category":"page"},{"location":"man/formatting/#Tables-with-Network-Dependent-Units","page":"Turning data into a CausalTable","title":"Tables with Network-Dependent Units","text":"","category":"section"},{"location":"man/formatting/","page":"Turning data into a CausalTable","title":"Turning data into a CausalTable","text":"The previous example assumes that each unit (row in the Table, in this case df), is \"causally independent\" of every other unit – that is, the treatment of one unit does not affect the response of any other unit. In some cases, however, we might work with data in which units may not be causally independent, but rather, in which one unit's variables could dependent on some summary function of its neighbors. ","category":"page"},{"location":"man/formatting/","page":"Turning data into a CausalTable","title":"Turning data into a CausalTable","text":"In this case, we can specify a graph argument to the CausalTable constructor, a Graph object from Graphs.jl which will be used to determine which units are neighbors of one another. We would also specify a summaries argument, a NamedTuple of NetworkSummary objects representing variables summarized over each unit's neighbors in the graph. More detail on the types of NetworkSummary that can be used in a dependent-data CausalTable can be found in Network Summaries","category":"page"},{"location":"man/formatting/","page":"Turning data into a CausalTable","title":"Turning data into a CausalTable","text":"Here's an example of how such a CausalTable might be constructed, using the Karate Club dataset. Treatment is defined as the number of friends a club member has, denoted by the summary function parameter summaries = (friends = Friends(),). ","category":"page"},{"location":"man/formatting/","page":"Turning data into a CausalTable","title":"Turning data into a CausalTable","text":"using CausalTables\nusing MLDatasets\nusing Graphs\n\n# Get a Table of Karate Club data from MLDatasets\ndata = KarateClub()\ntbl = data.graphs[1].node_data\n\n# Convert the karate club data into a Graphs.jl graph object\ng = SimpleGraphFromIterator([Edge(x...) for x in zip(data.graphs[1].edge_index...)])\n\n# Store the \"friends\" as an the adjacency matrix in a NamedTuple\n# Note that the input to arrays must be a NamedTuple, even if there is only one summary variable, \n# so the trailing comma is necessary.\nm = (F = adjacency_matrix(g),)\n\n# Construct a CausalTable with the adjacency matrix stored in `arrays` and a summary variable recording the number of friends\nctbl = CausalTable(tbl; treatment = :friends, response = :labels_clubs, arrays = m, summaries = (friends = Friends(:F),))\n\nnothing # hide","category":"page"},{"location":"#CausalTables.jl","page":"Home","title":"CausalTables.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A package for storing and simulating data for causal inference in Julia.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CausalTables.jl can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pkg> add CausalTables","category":"page"},{"location":"#Quick-Start","page":"Home","title":"Quick Start","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CausalTables.jl has three main functionalities:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Generating simulation data using a StructuralCausalModel\nComputing \"ground truth\" conditional distributions, means, and other functionals from a DataGeneratingProcess and a CausalTable\nWrapping an existing Table to make it a CausalTable for use by external packages.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The examples below illustrate each of these three functionalities.","category":"page"},{"location":"#Simulating-Data-from-a-DataGeneratingProcess","page":"Home","title":"Simulating Data from a DataGeneratingProcess","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To set up a statistical simulation using CausalTables.jl, we first define a StructuralCausalModel (SCM). This consists of two parts: a DataGeneratingProcess (DGP) that controls how the data is generated, and a list of variables to define the basic structure of the underlying causal diagram.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A DataGeneratingProcess can be constructed using the @dgp macro, which takes a sequence of conditional distributions of the form [variable name] ~ Distribution(args...) and returns a DataGeneratingProcess object. Then, one can construct an StructuralCausalModel by passing the DGP to its construct, along with labels of the treatment and response variables.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using CausalTables\nusing Random\nusing Distributions\n\ndgp = @dgp(\n        W ~ DiscreteUniform(1, 5),\n        X ~ (@. Normal(W, 1)),\n        Y ~ (@. Normal(X + 0.2 * W, 1))\n    )\n\nscm = StructuralCausalModel(\n    dgp;\n    treatment = :X,\n    response = :Y,\n    confounders = [:W]\n)","category":"page"},{"location":"","page":"Home","title":"Home","text":"One we've defined our list of distribution functions, we can generate data from the DGP using the rand function:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Random.seed!(1);\ndata = rand(scm, 5)","category":"page"},{"location":"","page":"Home","title":"Home","text":"For a more detailed guide of how to generate data please refer to Generating Data.","category":"page"},{"location":"#Computing-Ground-Truth-Conditional-Distributions","page":"Home","title":"Computing Ground Truth Conditional Distributions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Once we've defined a DGP and have some table of data with variables matching those of our DGP, we can compute the ground truth conditional distributions of any variable in a CausalTable (given a corresponding DGP) using the condensity function. This returns a Distribution object from the package Distributions.jl","category":"page"},{"location":"","page":"Home","title":"Home","text":"X_distribution = condensity(scm, data, :X)\n\n# output\n5-element Vector{Normal{Float64}}:\n Distributions.Normal{Float64}(μ=1.0, σ=1.0)\n Distributions.Normal{Float64}(μ=2.0, σ=1.0)\n Distributions.Normal{Float64}(μ=4.0, σ=1.0)\n Distributions.Normal{Float64}(μ=4.0, σ=1.0)\n Distributions.Normal{Float64}(μ=5.0, σ=1.0)","category":"page"},{"location":"","page":"Home","title":"Home","text":"For convenience, there also exists a conmean function that extracts the true conditional mean of a specific variable the CausalTable:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Y_mean = conmean(scm, data, :Y)\n\n# output\n5-element Vector{Float64}:\n 1.467564418628885\n 4.149933692528245\n 3.973979208080703\n 3.757247582108903\n 5.670866154143596","category":"page"},{"location":"","page":"Home","title":"Home","text":"For a more detailed guide of how to compute ground truth conditional distributions please refer to Computing Ground Truth Conditional Distributions.","category":"page"},{"location":"#Wrapping-an-existing-Table-as-a-CausalTable","page":"Home","title":"Wrapping an existing Table as a CausalTable","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you have a table of data that you would like to use with CausalTables.jl without defining a corresponding DataGeneratingProcess (i.e. to use with another package) you can wrap it as a CausalTable using the corresponding constructor:","category":"page"},{"location":"","page":"Home","title":"Home","text":"tbl = (W = rand(1:5, 10), X = randn(10), Y = randn(10))\nctbl = CausalTable(tbl; treatment = :X, response = :Y, confounders = [:W])","category":"page"},{"location":"","page":"Home","title":"Home","text":"For a more detailed guide of how to wrap an existing table as a CausalTable please refer to Turning Your Data Into a CausalTable.","category":"page"},{"location":"man/api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"man/api/","page":"API","title":"API","text":"The following documents public methods of CausalTables.jl.","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"Modules = [CausalTables]\nOrder   = [:function, :type]","category":"page"},{"location":"man/api/#Base.rand-Tuple{StructuralCausalModel, Int64}","page":"API","title":"Base.rand","text":"rand(scm::StructuralCausalModel, n::Int)\n\nGenerate random data from a Structural Causal Model (SCM) using the specified number of samples.\n\nArguments\n\nscm::StructuralCausalModel: The Structural Causal Model from which to generate data.\nn::Int: The number of samples to generate.\n\nReturns\n\nA CausalTable object containing the generated data.\n\n\n\n\n\n","category":"method"},{"location":"man/api/#CausalTables.condensity-Tuple{StructuralCausalModel, CausalTable, Symbol}","page":"API","title":"CausalTables.condensity","text":"condensity(scm::StructuralCausalModel, ct::CausalTable, var::Symbol)\n\nCompute the conditional density of a variable in a StructuralCausalModel given a CausalTable.\n\nArguments\n\nscm::StructuralCausalModel: The StructuralCausalModel representing the data generating process.\nct::CausalTable: The CausalTable containing the observed data.\nvar::Symbol: The variable for which to compute the conditional density.\n\nReturns\n\nThe conditional density of the variable var given the observed data.\n\n\n\n\n\n","category":"method"},{"location":"man/api/#CausalTables.conmean-Tuple{StructuralCausalModel, CausalTable, Symbol}","page":"API","title":"CausalTables.conmean","text":"conmean(dgp::DataGeneratingProcess, ct::CausalTable, var::Symbol)\n\nCompute the conditional mean of a variable in a CausalTable based on a DataGeneratingProcess.\n\nArguments\n\ndgp::DataGeneratingProcess: The DataGeneratingProcess object representing the data generating process.\nct::CausalTable: The CausalTable object representing the data.\nvar::Symbol: The variable for which to compute the conditional mean.\n\nReturns\n\nAn array of conditional means for the specified variable.\n\n\n\n\n\n","category":"method"},{"location":"man/api/#CausalTables.getscm-Tuple{CausalTable}","page":"API","title":"CausalTables.getscm","text":"getscm(o::CausalTable)\n\nGet the structural causal model (SCM) of a CausalTable object.\n\nThis function merges the column table of the CausalTable object with its arrays.\n\nArguments\n\no::CausalTable: The CausalTable object.\n\nReturns\n\nA merged table containing the column table and arrays of the CausalTable object.\n\n\n\n\n\n","category":"method"},{"location":"man/api/#CausalTables.replace-Tuple{CausalTable}","page":"API","title":"CausalTables.replace","text":"replace(o::CausalTable; kwargs...)\n\nReplace the fields of a CausalTable object with the provided keyword arguments.\n\nArguments\n\no::CausalTable: The CausalTable object to be replaced.\nkwargs...: Keyword arguments specifying the new values for the fields.\n\nReturns\n\nA new CausalTable object with the specified fields replaced.\n\n\n\n\n\n","category":"method"},{"location":"man/api/#CausalTables.summarize-Tuple{CausalTable}","page":"API","title":"CausalTables.summarize","text":"summarize(o::CausalTable; keep_original=true)\n\nSummarizes the data in a CausalTable object.\n\nArguments\n\no::CausalTable: The CausalTable object to be summarized.\nkeep_original::Bool: Whether to keep the original data in the resulting CausalTable. Default is true.\n\nReturns\n\nIf keep_original is true, a new CausalTable object with the original data merged with the summarized data.\nIf keep_original is false, a dictionary containing the summarized data.\n\n\n\n\n\n","category":"method"},{"location":"man/api/#Distributions.convolve-Union{Tuple{Vector{T}}, Tuple{T}} where T<:(Distributions.UnivariateDistribution)","page":"API","title":"Distributions.convolve","text":"Distributions.convolve(ds::Vector{T}) where {T <: UnivariateDistribution}\n\nOverload the convolve function to work on a vector of UnivariateDistribution.\n\nArguments\n\nds::Vector{T}: A vector of UnivariateDistribution objects.\n\nReturns\n\noutput: The result of convolving all the distributions in ds. If ds is empty, will return Binomial(0, 0.5) denoting a point mass at 0.\n\n\n\n\n\n","category":"method"},{"location":"man/api/#CausalTables.DataGeneratingProcess","page":"API","title":"CausalTables.DataGeneratingProcess","text":"mutable struct DataGeneratingProcess\n\nA struct representing a data generating process.\n\nFields\n\nnames: An array of symbols representing the names of the variables.\ntypes: An array of symbols representing the types of the variables.\nfuncs: An array of functions representing the generating functions for each variable.\n\n\n\n\n\n","category":"type"},{"location":"man/api/#CausalTables.Friends","page":"API","title":"CausalTables.Friends","text":"mutable struct Friends <: NetworkSummary\n\nA NetworkSummary counting the number of connected individuals in an adjacency matrix, also known as the number of \"friends\"\n\nFields\n\nmatrix::Symbol: The matrix representing the network.\n\n\n\n\n\n","category":"type"},{"location":"man/api/#CausalTables.NetworkSummary","page":"API","title":"CausalTables.NetworkSummary","text":"abstract type NetworkSummary\n\nAbstract type representing a summary of a network.\n\n\n\n\n\n","category":"type"},{"location":"man/api/#CausalTables.StructuralCausalModel","page":"API","title":"CausalTables.StructuralCausalModel","text":"struct StructuralCausalModel\n\nA struct representing a structural causal model (SCM).\n\nFields\n\ndgp::DataGeneratingProcess: The data generating process associated with the structural causal model.\ntreatment::Vector{Symbol}: The variables representing the treatment in the structural causal model.\nresponse::Vector{Symbol}: The variables representing the response in the structural causal model.\nconfounders::Vector{Symbol}: The variables representing the confounders in the structural causal model.\n\nConstructors\n\nStructuralCausalModel(dgp, treatment, response, confounders): Constructs a new StructuralCausalModel object.\n\nArguments\n\ndgp: The data generating process associated with the structural causal model.\ntreatment: The variables representing the treatment in the structural causal model.\nresponse: The variables representing the response in the structural causal model.\nconfounders: The variables representing the confounders in the structural causal model.\n\n\n\n\n\n","category":"type"},{"location":"man/api/#CausalTables.Sum","page":"API","title":"CausalTables.Sum","text":"Sum <: NetworkSummary\n\nA NetworkSummary which sums the values of the target variable for each unit connected in the adjacency matrix.\n\nFields\n\ntarget::Symbol: The target variable of the network.\nmatrix::Symbol: The matrix representation of the network.\n\n\n\n\n\n","category":"type"}]
}
