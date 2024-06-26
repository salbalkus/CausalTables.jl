# Turning Your Data Into a `CausalTable`

One of the main purposes of CausalTables.jl is to wrap a Table of data in Julia in order to provide it as input to some other causal inference package. Given a Table of some data, we can turn it into a `CausalTable` by specifying the treatment, response, and control variables. 

## Tables with Causally Independent Units

The code below demonstrates this on the Titanic dataset. This could be, for example, to use as input into some estimator of whether a passenger's sex caused them to survive the Titanic disaster, controlling for some baselineline confounders listed in `confounders`.

```@example titanic
using CausalTables
using MLDatasets: Titanic
using DataFrames

df = Titanic().dataframe

# Wrapping the dataset in a CausalTable
ctbl = CausalTable(df; treatment = :Sex, response = :Survived, confounders = [:Pclass, :Age, :SibSp])

nothing # hide
```

## Tables with Network-Dependent Units

The previous example assumes that each unit (row in the Table, in this case `df`), is "causally independent" of every other unit -- that is, the treatment of one unit does not affect the response of any other unit. In some cases, however, we might work with data in which units may *not* be causally independent, but rather, in which one unit's variables could dependent on some summary function of its neighbors. 

In this case, we can specify a `graph` argument to the `CausalTable` constructor, a Graph object from Graphs.jl which will be used to determine which units are neighbors of one another. We would also specify a `summaries` argument, a `NamedTuple` of `NetworkSummary` objects representing variables summarized over each unit's neighbors in the graph. More detail on the types of `NetworkSummary` that can be used in a dependent-data `CausalTable` can be found in [Network Summaries](network-summaries.md)

Here's an example of how such a `CausalTable` might be constructed, using the Karate Club dataset. Treatment is defined as the number of friends a club member has, denoted by the summary function parameter `summaries = (friends = Friends(),)`. 

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
m = (F = adjacency_matrix(g),)

# Construct a CausalTable with the adjacency matrix stored in `arrays` and a summary variable recording the number of friends
ctbl = CausalTable(tbl; treatment = :friends, response = :labels_clubs, arrays = m, summaries = (friends = Friends(:F),))

nothing # hide
```