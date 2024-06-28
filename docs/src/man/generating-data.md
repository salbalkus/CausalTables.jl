# Generating Data for Statistical Experiments

When evaluating a causal inference method, we often want to test it on data from a known causal model. CausalTables.jl allows us to define a `DataGeneratingProcess` (or DGP) to do just that. 

## Defining a DataGeneratingProcess

A data generating process describes a mechanism by which draws from random variables are simulated. It typically takes the form of a sequence of conditional distributions. CausalTables allows us to define a DGP as a DataGeneratingProcess object, which takes three arguments: the `names` of variables generated at each step, the `types` of these variables, and `funcs`, an array of functions of the form `(; O...) -> *some code`. 

Suppose, for example, that we wanted to simulate data from the following DGP:

```math
\begin{align*}
    W &\sim \text{DiscreteUniform}(1, 5) \\
    X &\sim \text{Normal}(W, 1) \\
    Y &\sim \text{Normal}(X + 0.2W, 1)
\end{align*}
```

where `X` is the treatment, `Y` is the response, and `W` is a confounding variable affecting both X and Y. A verbose and inconvenient (albeit correct) way to define this DGP would be as follows:

```jldoctest generation; output = false, filter = r"(?<=.{21}).*"s
using Distributions
using CausalTables

DataGeneratingProcess(
    [:W, :X, :Y],
    [:distribution, :distribution, :distribution],
    [
        (; O...) -> DiscreteUniform(1, 5), 
        (; O...) -> (@. Normal(O.W, 1)),
        (; O...) -> (@. Normal(O.X + 0.2 * O.W, 1))
    ]
)

# output
DataGeneratingProcess
```
where `; O...` syntax is a shorthand for a function that takes keyword arguments corresponding to the names of the variables in the DGP. 

However, a much more convenient way to define this DGP is using the `@dgp` macro, which takes a sequence of conditional distributions of the form `[variable name] ~ Distribution(args...)` and deterministic variable assignments of the form `[variable name] = f(...)` and automatically generates a valid DataGeneratingProcess. For example, the *easier* way to define the DGP above is as follows:

```jldoctest generation; output = false, filter = r"(?<=.{21}).*"s
using CausalTables
distributions = @dgp(
        W ~ DiscreteUniform(1, 5),
        X ~ (@. Normal(W, 1)),
        Y ~ (@. Normal(X + 0.2 * W, 1))
    )

# output
DataGeneratingProcess
```

Note that with the `@dgp` macro, any symbol (that is, any string of characters prefixed by a colon, as in `:W` or `:X`) is automatically replaced with the corresponding previously-defined variable in the process. For instance, in `Normal(:W, 1)`, the `:W` will be replaced automatically with the distribution we defined as `W` earlier in the sequence. 

## Defining a StructuralCausalModel

In CausalTables.jl, a StructuralCausalModel is a data generating process endowed with some causal interpretation. Constructing a StructuralCausalModel allows users to randomly draw a CausalTable with the necessary components from the DataGeneratingProcess they've defined. With the above DataGeneratingProcess in hand, we can define a `StructuralCausalModel` object like so -- treatment, response, and confounder variables in the causal model are specified as keyword arguments to the `DataGeneratingProcess` constructor:


```jldoctest generation; output = false, filter = r"(?<=.{21}).*"s
dgp = StructuralCausalModel(
    distributions;
    treatment = :X,
    response = :Y,
    confounders = [:W]
)

# output
StructuralCausalModel
```

## Networks of Causally-Connected Units

In some cases, we might work with data in which units may *not* be causally independent, but rather, in which one unit's variables could dependent on some summary function of its neighbors. Generating data from such a model can be done by adding lines of the form `Xs $ NetworkSummary` to the `@dgp` macro.

Here's an example of how such a `DataGeneratingProcess` might be constructed:

```jldoctest network; output = false, filter = r"(?<=.{21}).*"s
using Graphs
using CausalTables

dgp = @dgp(
        W ~ DiscreteUniform(1, 5),
        n = length(W),
        A = adjacency_matrix(erdos_renyi(n, 0.5)),
        Ws $ Sum(:W, :A),
        X ~ (@. Normal(Ws, 1)),
        Xs $ Sum(:X, :A),
        Y ~ (@. Normal(Xs + 0.2 * Ws, 1))
    )

scm = StructuralCausalModel(
    dgp;
    treatment = :X,
    response = :Y,
    confounders = [:W, :Ws]
)

# output
StructuralCausalModel
```



