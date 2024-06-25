macro scm(args...)
    # parse each line of the input into a vector of vectors
    parsed_result = [_parse_step(arg) for arg in args]

    # dump the vector of vectors into the DataGeneratingProcess constructor
    return DataGeneratingProcess(zip(parsed_result...))
end

mutable struct DataGeneratingProcess
    names::Vector{Union{Symbol, String}}
    funcs::Vector{Function}
end

DataGeneratingProcess(steps) = all(length(steps[1]) .== length.(arr[2:end])) ? DataGeneratingProcess(steps) : throw(ArgumentError("All step vectors must be the same length."))
length(x::DataGeneratingProcess) = length(x.names)

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

    elseif expr.args[1] == :(:=)
        # construct a function representing the step at the DGP
        func = eval(:((; O...) -> $(_replace_symbols_with_index(expr.args[3]))))
        
    elseif expr.args[1] == :(=)
        # construct a function representing the step at the DGP
        func = eval(:((; O...) -> $(expr.args[3])))

    else
        throw(ArgumentError("Invalid expression. Each line in the dgp macro must be of the form `var ~ distribution`, `var := code`, or `var = transformation``."))
    end

    return (name, func, can_sample, can_transform)
end

# Helper function to parse the dgp macro into anonymous functions
function _replace_symbols_with_index(expr)
    return postwalk(expr) do s
        typeof(s)==QuoteNode && return (:(O[$s]))
        s
    end
end

mutable struct StructuralCausalModel
    dgp::DataGeneratingProcess
    treatment::Vector{Symbol}
    response::Vector{Symbol}
    confounders::Vector{Symbol}

    function StructuralCausalModel(dgp, treatment, response, confounders)
        treatment, response, confounders = _process_causal_variable_names(treatment, response, confounders)
        new(dgp, treatment, response, confounders)
    end
end

length(scm::StructuralCausalModel) = length(scm.dgp)

# Function to draw a random CausalTable
function Base.rand(scm::StructuralCausalModel, n::Int)

    # Construct a Vector as subsequent inputs into each step function of the DGP
    result = Vector{Pair{Symbol, Any}}(undef, length(scm))
    tag = Vector{Symbol}(undef, length(scm))
    summaries = (;)

    # Iterate through each step of the SCM
    for i_step in 1:length(scm)

        # Draw from the result of the step function
        output_step = scm.funcs[i_step](output...)
        result[i_step] = scm.dgp.names[i_step] => _scm_draw(output_step, output, n)
        
        # Determine where the resulting output should be placed in the CausalTable
        if result[i_step] <: AbstractVector{<:Number}
            tag[i_step] = :data
        elseif result[i_step] <: AbstractArray{<:Number, 2}
            tag[i_step] = :arrays
        else
            tag[i_step] = :nothing
        end

        # If the initial output was a NetworkSummary, record this to include in the CausalTable
        if output_step <: NetworkSummary
            summaries = merge(summaries, NamedTuple{scm.name[i_step]}((output_step,)))
        end    
    end

    data = NamedTuple(result[tag .== :data])
    arrays = NamedTuple(result[tag .== :arrays])

    # Store the recorded draws in a CausalTable format
    return CausalTable(data, scm.treatment, scm.response, scm.confounders, arrays, summaries)
end

# Helper function for drawing a random CausalTable
_scm_draw(x::Distribution,            o::CausalTable, n::Int64) = rand(x, n)
_scm_draw(x::Vector{Distribution},    o::CausalTable, n::Int64) = length(x) == n ? rand.(x) : throw(ArgumentError("Length of vector of distributions in DataGeneratingProcess must be equal to n"))
_scm_draw(x::NetworkSummary,          o::CausalTable, n::Int64) = summarize(o, x)
_scm_draw(x::AbstractArray{<:Number}, o::CausalTable, n::Int64) = x


# Get the conditional density of a variable in the DGP
function condensity(scm::StructuralCausalModel, ct::CausalTable, var::Symbol)
    
end