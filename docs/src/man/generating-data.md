# Generating Data for Statistical Experiments

When evaluating a causal inference method, we often want to test it on data from a known causal model. CausalTables.jl allows us to define a `DataGeneratingProcess` (or DGP) to do just that. 

## Defining a DataGeneratingProcess

A data generating process describes a mechanism by which draws from random variables are simulated. It typically takes the form of a sequence of conditional distributions. CausalTables allows us to define a DGP as a DataGeneratingProcess object, which takes three arguments: the `names` of variables generated at each step, the `types` of these variables, and `funcs`, an array of functions of the form `O -> *some code`. 

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
    [
        O -> DiscreteUniform(1, 5), 
        O -> (@. Normal(O.W, 1)),
        O -> (@. Normal(O.X + 0.2 * O.W, 1))
    ]
)

# output
DataGeneratingProcess
```
where `O` is an object that stores the output of each previous function in the sequence as a field with a name corresponding to its order in the sequence (i.e. in this example, the first function's output is stored as `O.W`, the second function's output is stored as `O.X`, and so on).

However, a much more convenient way to define this DGP is using the `@dgp` macro, which takes a sequence of conditional distributions of the form `[variable name] ~ Distribution(args...)` and deterministic variable assignments of the form `[variable name] = f(...)` and automatically generates a valid DataGeneratingProcess. For example, the *easier* way to define the DGP above is as follows:

```jldoctest generation; output = false, filter = r"(?<=.{21}).*"s
distributions = @dgp(
        W ~ DiscreteUniform(1, 5),
        X ~ Normal.(W, 1),
        Y ~ (@. Normal(X + 0.2 * W, 1))
    )

# output
DataGeneratingProcess
```

Note that when using the `@dgp` macro, any symbol defined on the left side of an equation in the sequence can be used to pass in the output of a previous step on the right side. For example, in the above code, the symbol `W` is used to pass in the output of the first step to the second step. This works by metaprogramming which replaces `W` with `O.W` when the function is constructed by `@dgp`. 

We can also define steps other than distributions. There are four different types of "steps" that can be defined in a DGP sequence, each being constructed from a different "linking" symbol. Consider the following example, which uses all four types of steps:

```jldoctest generation; output = false, filter = r"(?<=.{21}).*"s
using Graphs
@dgp(
    W ~ Poisson(1),
    θ = exp.(W .+ 1),
    X ~ Normal.(θ, θ),
    G ≈ erdos_renyi(10, 0.5),
    M = Graphs.adjacency_matrix(ER),
    Xs $ Sum(:X, :G)
)

# output
DataGeneratingProcess
```

- Each `~` is used to denote a *Distribution* from Distributions.jl. These can both generate random data as well as admit expressions for the exact conditional distribution when calling functions like `condensity` (See [Computing ground truth conditional distributions](ground-truth.md)).
- The `=` symbol is used to denote *deterministic* functions of previous steps. They can be used to easily compute and reuse transformations of random variables. When a function like `condensity` is called on a `CausalTable`, each step will be recomputed to propagate any changes or interventions that may have been made, on the table.
- The `≈` symbol is used to denote *random* functions of previous steps that cannot necessarily be expressed as distributions -- for example, here we use `≈` to generate a random graph. When a function like `condensity` is called on a `CausalTable`, these steps will *not* be re-evaluated, so this symbol should not be used for functions depend on the values of previous steps.
- The `$` symbol is used to denote `NetworkSummary` functions. Similar to `=`, a NetworkSummary computes a deterministic transformation of previous steps, usually based on a random graph; the only difference is that when drawn from a StructuralCausalModel (see next section), the `NetworkSummary` will be stored in the CausalTable that is generated. See [Networks of Causally-Connected Units](#networks-of-causally-connected-units) or [Network summaries](network-summaries.md) for more details.

In this way, we can define virtually any DGP that can be expressed as a sequence of conditional distributions or transformations. For ease of use, one can still use the `O` object in the `@dgp` macro to pass in the output of all previous steps, which is especially useful for programmatically-defined DGPs. For example, the following code is equivalent to the above code:

```jldoctest generation; output = false, filter = r"(?<=.{21}).*"s
distributions = @dgp(
        W ~ DiscreteUniform(1, 5),
        X ~ Normal.(O[1], 1),
        Y ~ Normal.(hcat(values(O)...) * [1, 0.2], 1)
    )

# output
DataGeneratingProcess
```

In the first step, previous variables are accessed by index using `O[1]`, and in the third step, all previous variables are combined into a matrix by `hcat(values(O)...)`. Be careful when using these constructions, however, as they can make the code harder to read and understand. In some cases, it may be better to construct a `DataGeneratingProcess` manually using the constructor, for which several additional utilities are available. 

For instance, if one wanted to generate a large number of variables with the same distribution, one could use the `DataGeneratingProcess` constructor without specifying variable names, in which case names will be automatically generated:

```jldoctest generation; output = false, filter = r"(?<=.{75}).*"s
many_distributions = DataGeneratingProcess(
    [O -> Normal(0, 1) for _ in 1:100]
)

# output
DataGeneratingProcess([:X1, :X2, :X3, :X4, :X5, :X6, :X7, :X8, :X9, :X10  …
```

In addition, the `merge` function can be used to combine two separate DGP sequences into one:

```jldoctest generation; output = false, filter = r"(?<=.{21}).*"s
# Define a new distribution whose mean is the mean of previous draws
output_distribution = @dgp(
    Y ~ Normal.(reduce(+, values(O)) ./ n, 1)
)
# Merge our previous `many_distributions` with the new `output_distribution`
new_distributions = merge(many_distributions, output_distribution)

# output
DataGeneratingProcess
```

Finally, note that a DGP can depend on external variables. This is especially useful for running multiple simulations with different parameters, as one can define a function to generate DGPs from various sets of parameters:

```jldoctest generation; output = false, filter = r"(?<=.{21}).*"s
# Define a DGP that takes in parameters
dgp_family(a, b; σ2X = 1, σ2Y = 1) = @dgp(
        W ~ DiscreteUniform(a, b),
        X ~ Normal.(W, σ2X),
        Y ~ (@. Normal(X + 0.2 * W, σ2Y))
    )

# Create the same DGP but with different parameters
dgp_family(1, 5)
dgp_family(1, 10; σ2X = 2, σ2Y = 2)

# output
DataGeneratingProcess
```

Finally, if `dgp` denotes a `DataGeneratingProcess`, one can draw a sample path from it by calling `rand(dgp, n)` where `n` is the number of samples to draw. This will return a `NamedTuple` with the output of each step in the DGP. However, when running causal simulations, it is often more convenient to obtain a `CausalTable` object directly, which brings us to the next section: the `StructuralCausalModel`.

## Defining a StructuralCausalModel

In CausalTables.jl, a `StructuralCausalModel` is a data generating process endowed with some causal interpretation. Constructing a StructuralCausalModel allows users to randomly draw a CausalTable with the necessary components from the DataGeneratingProcess they've defined. With the previous DataGeneratingProcess in hand, we can define a `StructuralCausalModel` object like so -- treatment and response in the causal model are specified as keyword arguments to the `DataGeneratingProcess` constructor:


```jldoctest generation; output = false, filter = r"(?<=.{21}).*"s
scm = StructuralCausalModel(
    distributions;
    treatment = :X,
    response = :Y
)

# output
StructuralCausalModel
```

When a `StructuralCausalModel` is constructed with only treatment and response specified, all other variables are assumed to be confounders. However, one can also explicitly specify the causes of both treatment and response by passing them as a `NamedTuple` of lists to the `StructuralCausalModel` constructor:

```jldoctest generation; output = false, filter = r"(?<=.{21}).*"s
scm = StructuralCausalModel(
    distributions;
    treatment = :X,
    response = :Y,
    causes = (X = [:W], Y = [:X, :W])
)

# output
StructuralCausalModel
```

In the above, the keys of `causes` denote the variables whose causes are being specified, and the values are lists of variables that cause the key variable. In this case, the causes of the treatment `X` are specified as `[:W]`, and the causes of the response `Y` are specified as `[:X, :W]`, identical to how they are defined in a [CausalTable object](formatting.md). Just like for a `CausalTable`, while causes of other variables besides treatment and response can be specified, they are not necessary: only the causes of treatment and response are required as input. 

!!! note

    `causes` must be specified manually unless the user is assuming that all unlabeled variables cause both `treatment` and `outcome`. This is the default assumption of a `StructuralCausalModel`, but it may not not factually match the model encoded by the `DataGeneratingProcess`. This behavior is allowed for two reasons: (1) to permit a random draw of a `CausalTable` with an 'incorrect' causal model, which can be useful for benchmarking the robustness of different causal inference methods to model misspecification, and (2) to simulate causal models that implicitly condition on a particular set of variables by leaving them out of the `causes` argument. Otherwise, ensure that labels in `causes` do not contradict the data generating process! 

Finally, when setting up multiple simulations with similar DGPs and treatment/response labels, remember one can define a function to avoid repeating boilerplate code. Similar to how we defined a function earlier to generate multiple DGPs based on different sets of parameters, we can bundle everything together to create multiple SCMs:  

```jldoctest generation; output = false, filter = r"(?<=.{21}).*"s

scm_family(a, b; σ2X = 1, σ2Y = 1) = StructuralCausalModel(
    @dgp(
        W ~ DiscreteUniform(a, b),
        X ~ Normal.(W, σ2X),
        Y ~ (@. Normal(X + 0.2 * W, σ2Y))
    ); 
    treatment = :X,
    response = :Y
)

scm_family(1, 5)
scm_family(1, 10; σ2X = 2, σ2Y = 2)

# output
StructuralCausalModel
```

## Networks of Causally-Connected Units

In some cases, we might work with data in which units may *not* be causally independent, but rather, in which one unit's variables could dependent on some summary function of its neighbors. Generating data from such a model can be done by adding lines of the form `Xs $ NetworkSummary` to the `@dgp` macro.

Here's an example of how such a `StructuralCausalModel` might be constructed:

```jldoctest network; output = false, filter = r"(?<=.{21}).*"s
using Graphs
using CausalTables
using Distributions

dgp = @dgp(
        W ~ DiscreteUniform(1, 5),
        n = length(W),
        A = Graphs.adjacency_matrix(erdos_renyi(n, 0.5)),
        Ws $ Sum(:W, :A),
        X ~ (@. Normal(Ws, 1)),
        Xs $ Sum(:X, :A),
        Y ~ (@. Normal(Xs + 0.2 * Ws, 1))
    )

scm = StructuralCausalModel(
    dgp;
    treatment = :X,
    response = :Y
)

# output
StructuralCausalModel
```

## API

```@autodocs; canonical=false
Modules = [CausalTables]
Order   = [:type, :function]
Pages = ["data_generating_process.jl", "structural_causal_model.jl"]
```



