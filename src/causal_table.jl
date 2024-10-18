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
    # TODO: Undo temporary disabled line of code by checking for either 1 treatment/response or 2 treatments/response with 1 summarized
    #length(treatment) > 1 && throw(ArgumentError("Only univariate treatment is currently supported"))
    #length(response) > 1 && throw(ArgumentError("Only univariate response is currently supported"))

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
        #any(t ∉ names for t in treatment)   && throw(ArgumentError("Treatment variable(s) not found in data"))
        #any(r ∉ names for r in response)    && throw(ArgumentError("Response variable(s) not found in data"))
        #any(c ∉ names for c in confounders_replace_nothing) && throw(ArgumentError("Confounder variable(s) not found in data"))

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

_view_help(x::T, inds) where {T <: AbstractArray} = view(x, repeat([inds], ndims(x))...)
_view_help(x, inds) = x
_index_help(x::T, inds) where {T <: AbstractArray} = getindex(x, repeat([inds], ndims(x))...)
_index_help(x, inds) = x

function Tables.subset(o::CausalTable, inds; viewhint=nothing)
    viewhint = isnothing(viewhint) || viewhint
    
    data_subset = Tables.subset(o.data, inds; viewhint)

    if viewhint
        arrays_subset = map(x -> _view_help(x, inds), o.arrays)
    else
        arrays_subset = map(x -> _index_help(x, inds), o.arrays)
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
getscm(o::CausalTable) = merge(o.arrays, Tables.columntable(o.data)) # arrays must come first so that any summaries that are changed in the data are updated

Base.getindex(o::CausalTable, i::Int, j::Int) = Base.getindex(Tables.matrix(o.data), i, j)

function Base.show(io::IO, o::CausalTable)
    println(io, "CausalTable")
    PrettyTables.pretty_table(io, o, vcrop_mode=:middle, newline_at_end=false)
    println(io, "\nSummaries: $(o.summaries)")
    arrays_trunc = map(x -> typeof(x), o.arrays)
    println(io, "Arrays: $(arrays_trunc)")
end

# Functions to select causal variables from the data

"""
    select(o::CausalTable, symbols)

Selects specified columns from a `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which columns are to be selected.
- `symbols`: A list of symbols representing the columns to be selected.

# Returns
- A new `CausalTable` object with only the selected columns.

"""
select(o::CausalTable, symbols) = replace(o; data = o.data |> TableTransforms.Select(symbols...))

"""
    reject(o::CausalTable, symbols)

Removes the columns specified by `symbols` from the `CausalTable` object `o`.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which symbols will be rejected.
- `symbols`: A collection of symbols to be rejected from the `CausalTable`.

# Returns
A new `CausalTable` object with the specified symbols removed from its data.

"""
reject(o::CausalTable, symbols) = replace(o; data = o.data |> TableTransforms.Reject(symbols...))

"""
    treatment(o::CausalTable)

Selects the treatment column from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to select the treatment column.

# Returns
A new `CausalTable` containing only the treatment column
"""
treatment(o::CausalTable) = select(o, o.treatment)

"""
    confounders(o::CausalTable)

Selects and returns the confounders from a `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to select confounders.

# Returns
A new `CausalTable` containing only the confounders
"""
confounders(o::CausalTable) = select(o, o.confounders)

"""
    response(o::CausalTable)

Selects the response column from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to select the response column.

# Returns
A new `CausalTable` containing only the confounders
"""
response(o::CausalTable) = select(o, o.response)

"""
    treatmentparents(o::CausalTable)

Selects the confounders from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the parent variables of the treatment.

# Returns
A new `CausalTable` containing only the confounders
"""
treatmentparents(o::CausalTable) = reject(o, union(o.treatment, o.response))

"""
    responseparents(o::CausalTable)

Selects the treatment and confounders from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the parent variables of the response.

# Returns
A new `CausalTable` containing only the confounders and treatment
"""
responseparents(o::CausalTable) = reject(o, o.response)

# Other getters
"""
    data(o::CausalTable)

Retrieve the data stored in a `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` from which to retrieve the data.

# Returns
The data stored in the `CausalTable` object.
"""
data(o::CausalTable) = o.data

"""
    adjacency_matrix(o::CausalTable)

Generate the adjacency matrix induced by the `summaries` and `arrays` attributes of a `CausalTable` object. This matrix denotes which units are *causally dependent* upon one another: an entry of 1 in cell (i,j) indicates that some variable in unit i exhibits a causal relationship to some variable in unit j. 

# Arguments
- `o::CausalTable`: The `CausalTable` object for which the adjacency matrix is to be generated.

# Returns
A boolean matrix representing the adjacency relationships in the `CausalTable`.
"""
function adjacency_matrix(o::CausalTable)
    # Get the matrices used to summarize across observations in the table
    summary_matrix_names = unique([s.matrix for s in o.summaries if hasfield(typeof(s), :matrix)])
    if length(summary_matrix_names) > 0
        adj_matrices = values(o.arrays[summary_matrix_names])
        return(sum(adj_matrices) .!= 0.0)
    else
        return(LinearAlgebra.I(DataAPI.nrow(o)))
    end
end

"""
    dependency_matrix(o::CausalTable)

Generate the dependency matrix induced by the `summaries` and `arrays` attributes of a `CausalTable` object. This matrix stores which units are *statistically dependent* upon one another: an entry of 1 in cell (i,j) indicates that the data of unit i is correlated with the data in unit j. Two units are correlated if they either are causally dependent (neighbors in the adjacency matrix) or share a common cause (share a neighbor in the adjacency matrix).

# Arguments
- `o::CausalTable`: The `CausalTable` object for which the dependency matrix is to be generated.

# Returns
A boolean matrix representing the  relationships in the `CausalTable`.
"""
function dependency_matrix(o::CausalTable)
    # Get the matrices used to summarize across observations in the table
    summary_matrix_names = unique([s.matrix for s in o.summaries if hasfield(typeof(s), :matrix)])
    
    # Create a matrix where nonzero entries indicate that the two observations are dependent
    if length(summary_matrix_names) > 0

        # extract adjacency matrices from CausalTables
        adj_matrices = values(o.arrays[summary_matrix_names])

        # each unit to itself
        zero_hop = LinearAlgebra.I(DataAPI.nrow(o))

        # units that are neighbors
        one_hop = sum(adj_matrices)

        # units sharing a neighbor
        two_hop = sum(map(X -> X * X, adj_matrices))

        dependencies = (zero_hop + one_hop + two_hop) .> 0

    else
        dependencies = LinearAlgebra.I(DataAPI.nrow(o))
    end
    
    # map nonzero entires to 1
    return dependencies
end