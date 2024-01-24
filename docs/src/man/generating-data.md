# Generating Data for Statistical Experiments

When evaluating a causal inference method, we often want to test it on data from a known causal model. CausalTables.jl allows us to define a `DataGeneratingProcess` (or DGP) to do just that in two easy steps. First, we need to define a Vector of Pairs of the form `variable name => (; O...) -> Distribution(args...)`. The Distribution should take arguments corresponding to the names of the variables in the DGP. For example, suppose we wanted to define a statistical model of the form

\begin{align*}
    W &\sim \text{DiscreteUniform}(1, 5) \\
    X &\sim \text{Normal}(W, 1) \\
    Y &\sim \text{Normal}(X + 0.2W, 1)
\end{align*}

where `X` is the treatment, `Y` is the response, and `W` is a control variable. A verbose and inconvenient (albeit correct) way to define this DGP would be as follows:

```jldoctest generation
using Distributions

distributions = [
        :W => (; O...) -> DiscreteUniform(1, 5),
        :X => (; O...) -> (@. Normal(O[:L1], 1)),
        :Y => (; O...) -> (@. Normal(O[:A] + 0.2 * O[:L1], 1))
    ]
```

Note how Distributions can take previous variables as arguments by referencing them from the object `O`. The `; O...` syntax is a shorthand for a function that takes keyword arguments corresponding to the names of the variables in the DGP. 

However, a much more convenient way to define this DGP is using the `@dgp` macro, which takes a sequence of conditional distributions of the form `[variable name] ~ Distribution(args...)` and automatically generates a valid Vector of Pairs for a DataGeneratingProcess. For example, the *easier* way to define the DGP above is as follows:

```jldoctest generation
distributions = @dgp(
        W ~ DiscreteUniform(1, 5),
        X ~ (@. Normal(:W, 1)),
        Y ~ (@. Normal(:X + 0.2 * :W, 1))
    )

```

Note that with the `@dgp` macro, any symbol used in Distribution is automatically replaced with the corresponding previously-defined variable in the process. For instance, in `Normal(:W, 1)`, the `:W` will be replaced automatically with the distribution we defined as `W` earlier in the sequence. With this vector in hand, we can define a `DataGeneratingProcess` object like so -- treatment, response, and control variables in the causal model are specified as keyword arguments to the `DataGeneratingProcess` constructor:


```jldoctest generation
using CausalTables
dgp = DataGeneratingProcess(
    distributions;
    treatment = :X,
    response = :Y,
    controls = [:W]
)

```

## Networks of Causally-Connected Units


