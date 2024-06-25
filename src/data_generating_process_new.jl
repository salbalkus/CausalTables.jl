"""
    macro dgp(args...)

A macro to construct a DataGeneratingProcess from a sequence of distributions and transformations.

## Arguments
- `args...`: Variable number of arguments representing the steps of the data generating process.

## Returns
- `DataGeneratingProcess`: An instance of the `DataGeneratingProcess` type.

"""
macro dgp(args...)
    # parse each line of the input into a vector of vectors
    parsed_result = [_parse_step(arg) for arg in args]

    # dump the vector of vectors into the DataGeneratingProcess constructor
    return DataGeneratingProcess(zip(parsed_result...))
end

"""
    mutable struct DataGeneratingProcess

A struct representing a data generating process.

# Fields
- `names::Vector{Symbol}`: A vector of symbols representing the names of the variables.
- `types::Vector{Symbol}`: A vector of symbols representing the types of the variables.
- `funcs::Vector{Function}`: A vector of functions representing the generating functions for each variable.

"""
mutable struct DataGeneratingProcess
    names::Vector{Symbol}
    types::Vector{Symbol}
    funcs::Vector{Function}
end

DataGeneratingProcess(steps) = all(length(steps[1]) .== length.(arr[2:end])) ? DataGeneratingProcess(steps) : throw(ArgumentError("All step vectors must be the same length."))
Base.length(x::DataGeneratingProcess) = length(x.names)

# Helper function to parse each line in the dgp macro
function _parse_step(expr)

    ## 1) Parse the name
    name = expr.args[2]
    if !(name isa String || name isa Symbol)
        throw(ArgumentError("Invalid variable name. Variable name must be a string or symbol."))
    end

    ## 2) Parse the operation
    # `:=` means that the step is arbitrary code output that cannot be described in terms of a probability distribution
    # `=` means that the step is a transformation of a random variable
    # `~` indicates the step creates a distribution, which can either be sampled or used to compute a conditional density

    if expr.args[1] == :(~)
        # construct a function representing the step at the DGP
        func = eval(:((; O...) -> $(_replace_symbols_with_index(expr.args[3]))))
        type = :distribution

    elseif expr.args[1] == :(:=)
        # construct a function representing the step at the DGP
        func = eval(:((; O...) -> $(_replace_symbols_with_index(expr.args[3]))))
        type = :code
    elseif expr.args[1] == :(=)
        # construct a function representing the step at the DGP
        func = eval(:((; O...) -> $(expr.args[3])))
        type = :transformation
    else
        throw(ArgumentError("Invalid expression. Each line in the dgp macro must be of the form `var ~ distribution`, `var := code`, or `var = transformation``."))
    end

    return (name, func, type)
end

# Helper function to parse the dgp macro into anonymous functions
function _replace_symbols_with_index(expr)
    return postwalk(expr) do s
        typeof(s)==QuoteNode && return (:(O[$s]))
        s
    end
end