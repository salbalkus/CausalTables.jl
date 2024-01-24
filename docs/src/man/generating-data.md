# Generating Data for Statistical Experiments

When evaluating a causal inference method, we often want to test it on data from a known causal model. CausalTables.jl allows us to define a `DataGeneratingProcess` (or DGP) to do just that in two easy steps. First, we need to define a Vector of Pairs of the form `variable name => (; O...) -> Distribution(args...)`. The Distribution should take arguments corresponding to the names of the variables in the DGP. For example, suppose we wanted to define a statistical model of the form

\begin{align*}
    W &\sim \text{DiscreteUniform}(1, 5) \\
    X &\sim \text{Normal}(W, 1) \\
    Y &\sim \text{Normal}(X + 0.2W, 1)
\end{align*}

where `X` is the treatment, `Y` is the response, and `W` is a control variable. A verbose and inconvenient (albeit correct) way to define this DGP would be as follows:

```jldoctest generation; output = false, filter = r"(?<=.{16}).*"s
using Distributions

distributions = [
        :W => (; O...) -> DiscreteUniform(1, 5),
        :X => (; O...) -> (@. Normal(O[:W], 1)),
        :Y => (; O...) -> (@. Normal(O[:X] + 0.2 * O[:W], 1))
    ]

# output
3-element Vector
```

Note how Distributions can take previous variables as arguments by referencing them from the object `O`. The `; O...` syntax is a shorthand for a function that takes keyword arguments corresponding to the names of the variables in the DGP. 

However, a much more convenient way to define this DGP is using the `@dgp` macro, which takes a sequence of conditional distributions of the form `[variable name] ~ Distribution(args...)` and automatically generates a valid Vector of Pairs for a DataGeneratingProcess. For example, the *easier* way to define the DGP above is as follows:

```jldoctest generation; output = false, filter = r"(?<=.{16}).*"s
using CausalTables
distributions = @dgp(
        W ~ DiscreteUniform(1, 5),
        X ~ (@. Normal(:W, 1)),
        Y ~ (@. Normal(:X + 0.2 * :W, 1))
    )

# output
3-element Vector
```

Note that with the `@dgp` macro, any symbol used in Distribution is automatically replaced with the corresponding previously-defined variable in the process. For instance, in `Normal(:W, 1)`, the `:W` will be replaced automatically with the distribution we defined as `W` earlier in the sequence. With this vector in hand, we can define a `DataGeneratingProcess` object like so -- treatment, response, and control variables in the causal model are specified as keyword arguments to the `DataGeneratingProcess` constructor:


```jldoctest generation; output = false, filter = r"(?<=.{21}).*"s
dgp = DataGeneratingProcess(
    distributions;
    treatment = :X,
    response = :Y,
    controls = [:W]
)

# output
DataGeneratingProcess
```

## Networks of Causally-Connected Units

In some cases, we might work with data in which units may *not* be causally independent, but rather, in which one unit's variables could dependent on some summary function of its neighbors. Generating data from such a model can be done with two modifications:

1. Summary functions of neighbors are included in the `@dgp` macro via the form `varname = NetworkSummary(args...)`; see more information on avaible summaries in [Network Summaries](network-summaries.md)
2. The `DataGeneratingProcess` constructor is called with a function that takes in a sample size `n` and reutrns a `Graph` object from Graphs.jl, which is used to determine which units are neighbors of one another.

Here's an example of how such a `DataGeneratingProcess` might be constructed:

```jldoctest network; output = false, filter = r"(?<=.{21}).*"s
using Graphs
using CausalTables

distributions = @dgp(
        W ~ DiscreteUniform(1, 5),
        Ws = NeighborSum(:W),
        X ~ (@. Normal(:Ws, 1)),
        Xs = NeighborSum(:X),
        Y ~ (@. Normal(:Xs + 0.2 * :Ws, 1))
    )

dgp = DataGeneratingProcess(
    n -> erdos_renyi(n, 0.3),
    distributions;
    treatment = :Xs,
    response = :Y,
    controls = [:W, :Ws]
)

# output
DataGeneratingProcess
```



