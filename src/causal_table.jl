function _process_causal_variable_names(treatment, response, confounders)
    ## Process treatment and response variables into vectors, if they are not already vectors
    if typeof(treatment) <: Symbol
        treatment = [treatment]
    end

    if typeof(response) <: Symbol
        response = [response]
    end

    ## Error Handling ##

    # Assume only univariate treatment and response
    # TODO: Future code should allow for multivariate treatment and response
    length(treatment) != 1 && throw(ArgumentError("Only univariate treatment is currently supported"))
    length(response) != 1 && throw(ArgumentError("Only univariate response is currently supported"))

    # Ensure treatment, response, and confounders do not overlap
    name_occurrences = StatsBase.countmap(vcat(treatment, response, confounders))
    name_repeats = filter(((k, v),) -> v > 1, name_occurrences)
    length(name_repeats) > 0 && throw(ArgumentError("The following variable names are repeated across treatment, response, and confounder labels: $(keys(name_repeats))")) 
    
    # Return fully processed causal variable names
    return treatment, response, confounders
end

mutable struct CausalTable

    # data storage
    data::NamedTuple

    # labels
    treatment::Symbols
    response::Symbols
    confounders::Symbols

    # other
    arrays::NamedTuple
    summaries::NamedTuple

    function CausalTable(data, treatment, response, confounders, arrays, summaries)
        
        ## Process treatment and response variables into vectors, if they are not already vectors
        confounders_replace_nothing = isnothing(confounders) ? [] : confounders
        treatment, response, _ = _process_causal_variable_names(treatment, response, confounders_replace_nothing)

        # Ensure data input is a Table
        !Tables.istable(data) && throw(ArgumentError("`data` must be a Table. See https://tables.juliadata.org/ for more information."))        

        # Ensure treatment, response, and confounders are contained within the data
        names = (Tables.columnnames(Tables.columns(data))..., keys(summaries)...)
        any(t ∉ names for t in treatment)   && throw(ArgumentError("Treatment variable(s) not found in data"))
        any(r ∉ names for r in response)    && throw(ArgumentError("Response variable(s) not found in data"))
        any(c ∉ names for c in confounders_replace_nothing) && throw(ArgumentError("Confounder variable(s) not found in data"))

        # If confounders are Nothing, set them to be all columns (besides treatment and response) in the data by default
        if isnothing(confounders) 
            potential_confounders = union(Tables.columnnames(Tables.columns(data)), keys(summaries))
            not_confounders = vcat(treatment, response)
            confounders = setdiff(potential_confounders, not_confounders)
        end

        ## Construction ##

        # store the names of the input data columns
        names = Tables.columnnames(Tables.columns(data))  

        # store a matrix of data from the input table  
        data_table = Tables.columntable(data)
        
        # Construct a CausalTable with an underlying MatrixTable to store random vectors
        new(data_table, treatment, response, confounders, arrays, summaries)
    end
end

CausalTable(data, treatment, response; confounders = nothing, arrays = (;), summaries = (;)) = CausalTable(data, treatment, response, confounders, arrays, summaries)
CausalTable(data, treatment, response, confounders; arrays = (;), summaries = (;)) = CausalTable(data, treatment, response, confounders, arrays, summaries)
function CausalTable(data; treatment = nothing, response = nothing, confounders = nothing, arrays = (;), summaries = (;))
    isnothing(treatment) && throw(ArgumentError("Treatment variable must be defined"))
    isnothing(response) && throw(ArgumentError("Response variable must be defined"))
    CausalTable(data, treatment, response, confounders, arrays, summaries)
end

Tables.istable(::Type{CausalTable}) = true

### Column Interface ###
# Currently only allow column access from fixed data table (not network)

Tables.columnaccess(::Type{CausalTable}) = true
Tables.columns(o::CausalTable) = Tables.columns(o.data)

Tables.getcolumn(x::CausalTable, nm::Symbol) = Tables.getcolumn(Tables.columns(x.data), nm)
Tables.getcolumn(x::CausalTable, col::Int) = Tables.getcolumn(Tables.columns(x.data), col)
Tables.columnnames(x::CausalTable) = Tables.columnnames(Tables.columns(x.data))

# fixing StackOverflow error with column indexing methods via overloading
Tables.columnindex(x::CausalTable, nm::Symbol) = Tables.columnindex(x.data, nm)
Tables.columntype(x::CausalTable, nm::Symbol) = Tables.columntype(x.data, nm)

### Row Interface ###
# Currently only allow row access from fixed data table (not network)

Tables.rowaccess(::Type{CausalTable}) = true
rowaccess(::Type{<:CausalTable}) = true
rows(o::CausalTable) = Tables.rows(o.data)

### Other Tables Interface ###

Tables.schema(o::CausalTable) = Tables.schema(o.data)

# CausalTables do not permit materializers, because causal variable assignment is required via the constructor
#Tables.materializer(::Type{CausalTable})

function Tables.subset(o::CausalTable, inds; viewhint=nothing)
    viewhint = isnothing(viewhint) || viewhint
    
    data_subset = Tables.subset(o.data, inds; viewhint)

    if viewhint
        arrays_subset = map(x -> view(x, repeat([inds], ndims(x))...), o.arrays)
    else
        arrays_subset = map(x -> getindex(x, repeat([inds], ndims(x))...), o.arrays)
    end
    CausalTable(data_subset, o.treatment, o.response, o.confounders, arrays_subset, o.summaries)
end

DataAPI.nrow(o::CausalTable) = DataAPI.nrow(o.data)
DataAPI.ncol(o::CausalTable) = DataAPI.ncol(o.data)
### Causal-specific Features ###


"""
    replace(o::CausalTable; kwargs...)

Replace the fields of a `CausalTable` object with the provided keyword arguments.

# Arguments
- `o::CausalTable`: The `CausalTable` object to be replaced.
- `kwargs...`: Keyword arguments specifying the new values for the fields.

# Returns
A new `CausalTable` object with the specified fields replaced.

"""
replace(o::CausalTable; kwargs...) = CausalTable([field in keys(kwargs) ?  kwargs[field] : getfield(o, field) for field in fieldnames(typeof(o))]...)


"""
    getscm(o::CausalTable)

Get the structural causal model (SCM) of a `CausalTable` object.

This function merges the column table of the `CausalTable` object with its arrays.

# Arguments
- `o::CausalTable`: The `CausalTable` object.

# Returns
- A merged table containing the column table and arrays of the `CausalTable` object.

"""
getscm(o::CausalTable) = merge(Tables.columntable(o.data), o.arrays)

Base.getindex(o::CausalTable, i::Int, j::Int) = Base.getindex(Tables.matrix(o.data), i, j)

function Base.show(io::IO, o::CausalTable)
    println(io, "CausalTable")
    PrettyTables.pretty_table(io, o, vcrop_mode=:middle, newline_at_end=false)
    println(io, "\nSummaries: $(o.summaries)")
    arrays_trunc = map(x -> typeof(x), o.arrays)
    println(io, "Arrays: $(arrays_trunc)")
end

# Functions to get names of causal variables, including summarized versions
function _combine_summaries(o::CausalTable, symbols::AbstractArray{Symbol})
    sumnames = map(x -> gettarget(x), o.summaries)
    return vcat(symbols, [k for k in keys(sumnames)[values(sumnames) .∈ symbols]])
end

confoundernames(o::CausalTable) = _combine_summaries(o, o.confounders)
treatmentnames(o::CausalTable) = _combine_summaries(o, o.treatment)
responsenames(o::CausalTable) = _combine_summaries(o, o.response)

# Functions to select causal variables from the data
_summarize_select(o::CausalTable, symbols::AbstractArray{Symbol}) = summarize(o).data |> TableTransforms.Select(symbols...)

treatment(o::CausalTable) = _summarize_select(o, treatmentnames(o))
confounders(o::CausalTable) = _summarize_select(o, confoundernames(o))
response(o::CausalTable) = _summarize_select(o, responsenames(o))

# Functions to select the "previous" causal variables in the data generating process
treatmentparents(o::CausalTable) = summarize(o).data |> TableTransforms.Reject(treatmentnames(o)..., responsenames(o)...)
responseparents(o::CausalTable) = summarize(o).data |> TableTransforms.Reject(responsenames(o)...)


