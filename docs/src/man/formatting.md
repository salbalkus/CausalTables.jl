# Turning Your Data Into a `CausalTable`

In Julia, most datasets are stored in a Table: a data structure with a [Tables.jl](https://tables.juliadata.org/stable/)-compatible interface. One of the main purposes of CausalTables.jl is to wrap a Table of data in Julia in order to provide it as input to some other causal inference package. Given a Table of some data, we can turn it into a `CausalTable` by specifying the treatment, response, and control variables. 

## Tables with Causally Independent Units

The code below provides an example of how to wrap the Boston Housing dataset as a `CausalTable` to answer causal questions of the form "How would changing nitrous oxide air pollution (`NOX`) within Boston-area towns affect median home value (`MEDV`)?" Any dataset in a [Tables.jl](https://tables.juliadata.org/stable/)-compliant format can be wrapped as a `CausalTable`. In this example, we turn a `DataFrame` from [DataFrames.jl](https://dataframes.juliadata.org/stable/) into a `CausalTable` object.

```@example bostonhousing
using CausalTables
using MLDatasets: BostonHousing
using DataFrames

# get data in a Tables.jl-compliant format
tbl = BostonHousing().dataframe

# Wrapping the dataset in a CausalTable
ctbl = CausalTable(tbl; treatment = :NOX, response = :MEDV, confounders = [:CRIM, :ZN, :INDUS, :CHAS, :B, :DIS, :LSTAT])

nothing # hide
```

After wrapping a dataset in a `CausalTable` object, the [Tables.jl](https://tables.juliadata.org/stable/) is available to call on the `CausalTable` as well. Below, we demonstrate a few of these functions, as well as additional utility functions for causal inference tasks made available by CausalTables.jl.

```@example bostonhousing
using Tables

# Examples of using the Tables.jl interface
Tables.getcolumn(ctbl, :NOX) # extract specific column
Tables.subset(ctbl, 1:5)     # exact specific rows
Tables.columnnames(ctbl)     # obtain all column names

# Additional utility functions for CausalTables
treatment(ctbl)              # get CausalTable of treatment variables
response(ctbl)               # get CausalTable of response variables
confounders(ctbl)            # get CausalTable of confounders
responseparents(ctbl)        # get CausalTable of treatment and confounders
data(ctbl)                   # get underlying wrapped dataset

# replace one or more attributes of the CausalTable
CausalTables.replace(ctbl; response = :CRIM, confounders = [:MEDV, :ZN, :INDUS, :CHAS, :B, :DIS, :LSTAT]) 

nothing # hide
```

## Tables with Network-Dependent Units

The previous example assumes that each unit (row in the Table, in this case `tbl`), is "causally independent" of every other unit -- that is, the treatment of one unit does not affect the response of any other unit. This is a component of the "stable unit treatment value assumption" (SUTVA) often used in causal inference. In some cases, however, we might work with data in which units may *not* be causally independent, but rather, in which one unit's variables depend on some summary function of its neighbors

In this case, one must instead perform causal inference on the summary functions of each unit's neighbors ([Aronow and Samii, 2017](https://doi.org/10.1214/16-AOAS1005)). To do this, each `CausalTable` has two relevant arguments that can be used to correct SUTVA violations. The `arrays` argument is a `NamedTuple` that can store adjacency matrices and other miscellaneous parameters that denote the causal relationships between variables. The `summaries` argument is a tuple of `NetworkSummary` objects that can be used to summarize the network relationships between units by referencing variables in either the underlying data or the `arrays` argument of `CausalTable` (or both). 

The code below provides an example of how such a `CausalTable` might be constructed to consider a summary function treatment in the case of causally-dependent units, using the Karate Club dataset. In this example, treatment is defined as the number of friends a club member has, denoted by the summary function parameter `summaries = (friends = Friends(:F),)`. Hence, this answers the causal question "how would changing a subject's number of friends (`friends`) affect which club they are likely to join (`labels_clubs`)?" 

We store the network relationships between units as an adjacency matrix `F` by assigning it to the `arrays` parameters. This allows the `Friends(:F)` summary function to access it when calling `summarize(ctbl)`. More detail on the types of `NetworkSummary` that can be used in a dependent-data `CausalTable` can be found in [Network Summaries](network-summaries.md)

```@example karateclub
using CausalTables
using MLDatasets
using Graphs

# Get a Table of Karate Club data from MLDatasets
data = KarateClub()
tbl = data.graphs[1].node_data

# Convert the karate club data into a Graphs.jl graph object
g = SimpleGraphFromIterator([Edge(x...) for x in zip(data.graphs[1].edge_index...)])

# Store the "friends" as an the adjacency matrix in a NamedTuple
# Note that the input to arrays must be a NamedTuple, even if there is only one summary variable, 
# so the trailing comma is necessary.
m = (F = Graphs.adjacency_matrix(g),)

# Construct a CausalTable with the adjacency matrix stored in `arrays` and a summary variable recording the number of friends
ctbl = CausalTable(tbl; treatment = :friends, response = :labels_clubs, arrays = m, summaries = (friends = Friends(:F),))

nothing # hide
```

One can then call the function `summarize(ctbl)` to compute the values of the summary function on the causal table. 

Based on these summaries, it is also possible to extract two matrices from the `CausalTable` object: the `adjacency_matrix` and the `dependency_matrix`. The `adjacency_matrix` denotes which units are *causally dependent* upon one another: an entry of 1 in cell (i,j) indicates that some variable in unit i exhibits a causal relationship to some variable in unit j. The `dependency_matrix` stores which units are *statistically dependent* upon one another: an entry of 1 in cell (i,j) indicates that the data of unit i is correlated with the data in unit j. Two units are correlated if they either are causally dependent (neighbors in the adjacency matrix) or share a common neighbor in the adjacency matrix.

```@example karateclub
CausalTables.adjacency_matrix(ctbl) # get adjacency matrix
CausalTables.dependency_matrix(ctbl) # get dependency matrix

nothing # hide
```

## API for Working with a `CausalTable`
```@autodocs; canonical=false
Modules = [CausalTables]
Order   = [:type, :function]
Pages = ["causal_table.jl"]
```