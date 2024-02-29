# Turning Your Data Into a `CausalTable`

One of the main purposes of CausalTables.jl is to wrap a Table of data in Julia in order to provide it as input to some other causal inference package. Given a Table of some data, we can turn it into a `CausalTable` by specifying the treatment, response, and control variables. 

## Tables with Causally Independent Units

The code below demonstrates this on the Titanic dataset. This could be, for example, to use as input into some estimator of whether a passenger's sex caused them to survive the Titanic disaster, controlling for some baselineline covariates listed in `controls`.

```@example titanic
using CausalTables
using MLDatasets: Titanic
using DataFrames

df = Titanic().dataframe

# Wrapping the dataset in a CausalTable
ctbl = CausalTable(df; treatment = :Sex, response = :Survived, controls = [:Pclass, :Age, :SibSp])

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

# Note that the input to summaries must be a NamedTuple, even if there is only one summary variable, so the trailing comma is necessary.
ctbl = CausalTable(tbl; graph = g, treatment = :friends, response = :labels_clubs, summaries = (friends = Friends(),))

nothing # hide
```

Be warned: if you try to call `gettreatment` on a `CausalTable` that has not been summarized, you will get an error:

```@example karateclub
try  #hide
gettreatment(ctbl)
catch err; showerror(stderr, err); end  #hide
```

If you wish to extract the treatment variable, you will first need to call `summarize` on the CausalTable object, which computes the summary variables over the network. Then, calling `gettreatment` will yield the summarized treatment variable, like so:

```@example karateclub
ctbl_summarized = summarize(ctbl)
gettreatment(ctbl_summarized)
```

The response and controls can also be extracted using `getresponse` and `getcontrols`, respectively. 