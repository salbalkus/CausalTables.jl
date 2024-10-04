
"""
    macro dgp(args...)

A macro to construct a DataGeneratingProcess (DGP) from a sequence of distributions and transformations. A data generating process is a sequence of steps that generates a dataset. It can be used to encode the causal structure of a given statistical problem; for instance, if \$Y=f(X)\$ where \$f\$ is some function of \$X\$, then it can be said that \$X\$ *causes* \$Y\$. 

Each line in the `dgp` macro is treated as a discrete step of the DGP. In `CausalTables.jl`, a DGP object is required as an input to a StructuralCausalModel object. These generally have two uses: (a) randomly generating a dataset with a particular causal structure, and (b) computing the ground truth value of functions of the data. 

Each line generates one variable in the dataset using an assignment operator. When `rand` is called on a StructuralCausalModel that contains a DGP constructed via this macro, each line is computed in sequence depending on the assignment operator. Three assignment operators are available:

1. The standard `=` symbol, which evaluates the value of the right-hand side directly. This often denotes a **fixed** variable; i.e. `X = 0.5`. 
2. The `~` symbol, which constructs a Distribution object. This is used to denote a **random** variable; i.e. `X ~ Normal()`.
3. The `\$` symbol, which denotes a **summary** of random variables in previous steps; i.e. `X \$ Friends(:A)`.

Of course, one can always randomly generate random variables by calling a function on the right-hand side of `=`. The `~` operator serves two purposes. First, it allows a DGP to be expressed more concisely, by allowing just the distribution to be specified instead of needing to call `rand` at each step. Second, it allows the analytic computation of closed-form conditional densities, which can then be used to compute the ground truth value of statistical functionals that depend on the DGP. 

Lines with `\$` tell the `rand` function to generate the random variable at that step by summarizing the vector of previous vectors according to the NetworkSummary object on the right-hand side. Again, while this could be performed with `=`, the `\$` operator allows the closed-form conditional density of the summary function to be computed for certain summaries and random variables, such as sums of normal distributions. 

## Example

## Arguments
- `args...`: Variable number of arguments representing the steps of the data generating process.

## Returns
- `DataGeneratingProcess`: An instance of the `DataGeneratingProcess` type.

"""
macro dgp(args...)
    names = [_parse_name(arg) for arg in args]
    # parse each line of the input into a vector of vectors
    parsed_result = [_parse_step(arg, names[1:i]) for (i, arg) in enumerate(args)]

    # dump the vector of vectors into the DataGeneratingProcess constructor
    types, funcs  = collect(map(collect, zip(parsed_result...)))
    funcs = :([$(funcs...)])
    return :(DataGeneratingProcess($(names), $(types), $(esc(funcs))))
end

"""
    mutable struct DataGeneratingProcess

A struct representing a data generating process.

# Fields
- `names`: An array of symbols representing the names of the variables.
- `types`: An array of symbols representing the types of the variables.
- `funcs`: An array of functions representing the generating functions for each variable.

"""
mutable struct DataGeneratingProcess
    names::Symbols
    types::Symbols
    funcs::AbstractArray{Function, 1}
end

DataGeneratingProcess(steps) = all(length(steps[1]) .== length.(arr[2:end])) ? DataGeneratingProcess(steps) : throw(ArgumentError("All step vectors must be the same length."))
Base.length(x::DataGeneratingProcess) = length(x.names)

function _parse_name(expr)
    # Get the first value in the expression
    name = expr.args[length(expr.args)-1]

    # Test if the name is a string or symbol
    if !(name isa String || name isa Symbol)
        throw(ArgumentError("Invalid variable name. Variable name must be a string or symbol."))
    end

    return name
end

# Helper function to parse each line in the dgp macro
function _parse_step(expr, names)

    # `:=` means that the step is arbitrary code output that cannot be described in terms of a probability distribution
    # `=` means that the step is a transformation of a random variable
    # `~` indicates the step creates a distribution, which can either be sampled or used to compute a conditional density

    if length(expr.args) == 3 && expr.args[1] == :(~)
        # construct a function representing the step at the DGP
        func = :(O -> $(_replace_symbols_with_index(expr.args[3], names)))
        type = :distribution
    elseif length(expr.args) == 2
        # construct a function representing the step at the DGP
        func = :(O -> $(_replace_symbols_with_index(expr.args[2], names)))
        type = :code
    elseif length(expr.args) == 3 && expr.args[1] == :($)
        # construct a function representing the step at the DGP
        func = :(O -> $(expr.args[3]))
        type = :transformation
    else
        throw(ArgumentError("Invalid expression. Each line in the dgp macro must be of the form `var ~ distribution`, `var = code`, or `var \$ transformation`"))
    end

    return (type, func)
end

# Helper function to parse the dgp macro into anonymous functions
function _replace_symbols_with_index(expr, names)
    return postwalk(expr) do s
        if typeof(s)==Symbol
            if s in names
                return (:(O.$s))
            end
        end
        s
    end
end