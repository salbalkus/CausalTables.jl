function _process_causal_variable_names(treatment, response, causes)

    # Wrap symbols in lists if they haven't been already
    if treatment isa Symbol
        treatment = [treatment]
    end

    if response isa Symbol
        response = [response]
    end

    # Ensure treatment and response do not overlap
    name_occurrences = StatsBase.countmap(vcat(treatment, response))
    name_repeats = filter(((k, v),) -> v > 1, name_occurrences)
    length(name_repeats) > 0 && throw(ArgumentError("The following variable names are repeated across treatment and response lists: $(keys(name_repeats))")) 
    
    if !isnothing(causes)
        # Check that `causes` is acyclic
        _check_dag(causes) && throw(ArgumentError("`causes` contains a cycle, but causal relationships must form a directed acyclic graph (DAG), meaning no cycles are allowed."))

        # Check that `causes` includes all treatment and response variables
        !all(map(x -> x ∈ keys(causes), vcat(treatment, response))) && throw(ArgumentError("`causes` must contain, at a minimum, a key for each Symbol contained in `treatment` or `response`. If a given treatment has no causes, set it equal to an empty Vector (i.e. `A = []`). If this error arose from using the `replace` function to change the `treatment` or `response` of a `CausalTable`, be sure to also update the `causes` attribute to reflect the new changes."))

        # Ensure no responses cause any treatments in the DAG
        treatment_causes = map(t -> causes[t], treatment)
        response_is_a_cause = any.(map(r -> r .∈ treatment_causes, response))
        any(response_is_a_cause) && throw(ArgumentError("`causes` denotes one or more response variables causing a treatment variable, which is not allowed."))
    end

    # Return fully processed causal variable names
    return treatment, response, causes
end

# Function that outputs true if the input is not a directed acyclic graph
# Uses topological sorting to check
function _check_dag(causes_original)
    causes = deepcopy(causes_original)

    S = setdiff(keys(causes), vcat(values(causes)...))
    L = Set{Symbol}()

    while(length(S) > 0)
        s = pop!(S)
        push!(L, s)
        if s ∈ keys(causes)
            while length(causes[s]) > 0
                e = pop!(causes[s])
                if(e ∈ keys(causes))
                    push!(S, e)
                end
            end
        end
    end
    return any(length.(values(causes)) .!= 0)
end

# Default assumption: if not specified, set anything not labeled as a cause of all treatments and all responses
function set_unlabeled_causes(data, treatment, response)
    u = vcat(treatment, response)
    confounders = setdiff(Tables.columnnames(data), u)
    causes = Dict()

    for t in treatment
        causes[t] = confounders
    end

    confounders_and_treatment = vcat(confounders, treatment)
    for r in response
        causes[r] = confounders_and_treatment
    end
    return (;causes...) # Unpack dictionary into NamedTuple
end

mutable struct CausalTable
    # data storage
    data::NamedTuple

    # labels
    treatment::Symbols
    response::Symbols
    causes::NamedTuple

    # other
    arrays::NamedTuple
    summaries::NamedTuple

    function CausalTable(data, treatment, response, causes, arrays, summaries)

        # Ensure data input is a Table
        !Tables.istable(data) && throw(ArgumentError("`data` must be a Table. See https://tables.juliadata.org/ for more information."))

        # Store a matrix of data from the input table  
        data_table = Tables.columntable(data)

        ## Process treatment and response variables into vectors
        treatment, response, _ = _process_causal_variable_names(treatment, response, causes)        

        # Decide what to do when no causes are provided
        if(isnothing(causes))
            causes = set_unlabeled_causes(data, treatment, response)
        end
        
        # Construct a CausalTable with an underlying MatrixTable to store random vectors
        new(data_table, treatment, response, causes, arrays, summaries)
    end
end

CausalTable(data, treatment, response, causes; arrays = (;), summaries = (;)) = CausalTable(data, treatment, response, causes, arrays, summaries)
CausalTable(data, treatment, response; causes = nothing, arrays = (;), summaries = (;)) = CausalTable(data, treatment, response, causes, arrays, summaries)

function CausalTable(data; treatment = nothing, response = nothing, causes = nothing, arrays = (;), summaries = (;))
    isnothing(treatment) && throw(ArgumentError("Treatment variable must be defined"))
    isnothing(response) && throw(ArgumentError("Response variable must be defined"))
    CausalTable(data, treatment, response, causes, arrays, summaries)
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
    CausalTable(data_subset, o.treatment, o.response, o.causes, arrays_subset, o.summaries)
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

# Functions to select variables from the data

"""
    select(o::CausalTable, symbols)

Selects specified columns from a `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which columns are to be selected.
- `symbols`: A list of symbols representing the columns to be selected.

# Returns
- A new `CausalTable` object with only the selected columns.

"""
select(o::CausalTable, symbols::Symbol) = replace(o; data = o.data |> TableTransforms.Select(symbols))
select(o::CausalTable, symbols) = replace(o; data = (length(symbols) == 0 ? (;) : o.data |> TableTransforms.Select(symbols...)))

"""
    reject(o::CausalTable, symbols)

Removes the columns specified by `symbols` from the `CausalTable` object `o`.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which symbols will be rejected.
- `symbols`: A collection of symbols to be rejected from the `CausalTable`.

# Returns
A new `CausalTable` object with the specified symbols removed from its data.

"""
reject(o::CausalTable, symbols::Symbol) = replace(o; data = o.data |> TableTransforms.Reject(symbols))
reject(o::CausalTable, symbols) = replace(o; data = o.data |> TableTransforms.Reject(symbols...))


"""
    treatment(o::CausalTable)

Selects the treatment column(s) from the given `CausalTable` object.
treatment
# Arguments
- `o::CausalTable`: The `CausalTable` object from which to select the treatment column(s).

# Returns
A new `CausalTable` containing only the treatment column(s)
"""
treatment(o::CausalTable) = select(o, o.treatment)

"""
    treatmentmatrix(o::CausalTable)

Outputs the treatment column(s) from the given `CausalTable` object as a matrix.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to select the treatment column(s).

# Returns
A matrix containing only the treatment column(s)
"""
treatmentmatrix(o::CausalTable) = Tables.matrix(treatment(o))

"""
    response(o::CausalTable)

Selects the response column(s) from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to select the response column(s).

# Returns
A new `CausalTable` containing only the response column(s).
"""
response(o::CausalTable) = select(o, o.response)

"""
    responsematrix(o::CausalTable)

Outputs the response column(s) from the given `CausalTable` object as a matrix.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to select the response column(s).

# Returns
A matrix containing only the response column(s)
"""
responsematrix(o::CausalTable) = Tables.matrix(response(o))

"""
    parents(o::CausalTable, symbol)

Selects the variables that precede `symbol` causally from the CausalTable `o`, based on the `causes` attribute. Note that if `symbol` is not contained within `o.causes`, this function will output an empty `CausalTable`.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the parent variables of `symbol`.
- `symbol`: The variable for which to extract the parent variables.

# Returns
A new `CausalTable` containing only the parents of `symbol`
"""
function parents(o::CausalTable, symbol)
    if symbol in keys(o.causes)
        select(o, o.causes[symbol])
    else
        replace(o; data = (;))
    end
end

"""
    treatmentparents(o::CausalTable)

Selects the parents of each treatment variable from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the parent variables of each treatment.
- `collape_parents::Bool`: Optional parameter, whether to collapse the output to a single `CausalTable` object if there is either only one treatment or all treatments have the same parents. Defaults to `true`.

# Returns
A new `CausalTable` containing only the causes of the treatment (if a single treatment, or all treatments share the same set of causes); otherwise, a Vector of CausalTable objects containing the causes of each treatment.
"""
function treatmentparents(o::CausalTable; collapse_parents = true)
    # Extract the causes of each treatment variable as vector of vectors
    parent_names = [o.causes[k] for k in keys(o.causes) if k ∈ o.treatment]
    # When possible, only select a single `CausalTable` representing the parents of all variables
    if(collapse_parents && (length(o.treatment) == 1 || all(x -> x == parent_names[1], parent_names)))
        return(select(o, o.causes[o.treatment[1]]))
    # Otherwise, return a list of `CausalTables` representing the parents of each treatment
    else
        return([select(o, pn) for pn in parent_names])
    end
end


"""
    responseparents(o::CausalTable)

Selects the parents of each response variable from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the parent variables of each response.
- `collape_parents::Bool`: Optional parameter, whether to collapse the output to a single `CausalTable` object if there is either only one response or all response have the same parents. Defaults to `true`.

# Returns
A new `CausalTable` containing only the causes of the responses (if a single response, or all responses share the same set of causes); otherwise, a Vector of CausalTable objects containing the causes of each response.
"""
function responseparents(o::CausalTable; collapse_parents = true)
    # Extract the causes of each treatment variable as vector of vectors
    parent_names = [o.causes[k] for k in keys(o.causes) if k ∈ o.response]
    # When possible, only select a single `CausalTable` representing the parents of all variables
    if(collapse_parents && (length(o.response) == 1 || all(x -> x == parent_names[1], parent_names)))
        return(select(o, o.causes[o.response[1]]))
    # Otherwise, return a list of `CausalTables` representing the parents of each treatment
    else
        return([select(o, pn) for pn in parent_names])
    end
end

### Matrix utilities for the following

# Iterator over matrix of CausalTables to turn them into matrices
# If a CausalTable has empty data, return an empty matrix
matrix(tbl::Matrix{CausalTable}) = map(x -> length(x.data) == 0 ? [;] : Tables.matrix(x), tbl)
matrix(tbl::CausalTable) = Tables.matrix(tbl)

function select_over_matrix(o::CausalTable, varnames; collapse_parents = true) 
    # When possible, only select a single `CausalTable` representing the selection of all variables
    if(collapse_parents && (length(o.response) + length(o.treatment) == 1 || all(x -> x == varnames[1,1], vec(varnames))))
        return(select(o, varnames[1,1]))
    # Otherwise, return a matrix of `CausalTables` representing the selection of each treatment-response pair
    else
        return([select(o, cn) for cn in varnames])
    end
end

### Confounders ###

"""
    confoundernames(o::CausalTable)

Outputs the confounder names of each response-treatment pair from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the confounder names of each treatment-response pair.

# Returns
A matrix of Vectors containing the confounder names of each treatment-response pair.
"""
function confoundernames(o::CausalTable)
    # Extract the causes of each treatment variable as vector of vectors
    treatment_parent_names = [o.causes[k] for k in keys(o.causes) if k ∈ o.treatment]
    response_parent_names = [o.causes[k] for k in keys(o.causes) if k ∈ o.response]
    return(hcat([[intersect(tpn, rpn) for tpn in treatment_parent_names] for rpn in response_parent_names]...))
end

"""
    confoundernames(o::CausalTable, x::Symbol, y::Symbol)

Outputs the names of the confounders of the causal relationship between `x` and `y` from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the confounder names.
- `x::Symbol`, `y::Symbol`: The two variables whose confounders should be selected.

# Returns
A Vector of Symbols containing the names of the confounders between x and y.
"""
confoundernames(o::CausalTable, x::Symbol, y::Symbol) = intersect(o.causes[x], o.causes[y])

"""
    confounders(o::CausalTable; collapse_parents = true)

Selects the confounders of each response-treatment pair from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the confounder variables of each treatment-response pair.
- `collape_parents::Bool`: Optional parameter, whether to collapse the output to a single `CausalTable` object if there is either only one treatment-response pair or all pair share the same set of confounders. Defaults to `true`.

# Returns
A new `CausalTable` containing only the confounders (if a single response, or all responses share the same set of causes); otherwise, a Matrix of CausalTable objects containing the confounders of each treatment-response pair, where rows represent responses and columns represent treatments.
"""
confounders(o::CausalTable; collapse_parents = true) = select_over_matrix(o, confoundernames(o); collapse_parents = collapse_parents)

"""
    confounders(o::CausalTable, x::Symbol, y::Symbol)

Selects the common causes for a specific pair of variables (x,y) from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the confounders.
- `x::Symbol`, `y::Symbol`: The two variables whose confounders should be selected.

# Returns
A new `CausalTable` containing only the confounders of both x and y.
"""
confounders(o::CausalTable, x::Symbol, y::Symbol) = select(o, confoundernames(o, x, y))

"""
    confoundersmatrix(o::CausalTable; collapse_parents = true)

Outputs the treatment-variable confounders from the given `CausalTable` object as a matrix (or matrix of matrices, if multiple treatment-response pairs are present).

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the confounders of each treatment-response pair.
- `collape_parents::Bool`: Optional parameter, whether to collapse the output to a single `Matrix` object if there is either only one treatment-response pair or all pair share the same set of confounders. Defaults to `true`.

# Returns
A matrix containing only the confounders.
"""
confoundersmatrix(o::CausalTable; collapse_parents = true) = matrix(confounders(o; collapse_parents = collapse_parents))

### Mediation ###

"""
    mediatornames(o::CausalTable)

Outputs the mediator names of each response-treatment pair from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the mediator names of each treatment-response pair.

# Returns
A matrix of Vectors containing the mediator names of each treatment-response pair.
"""
function mediatornames(o::CausalTable)
    # Extract the mediator names: variables that cause response but are caused by treatment
    caused_by_treatment = [[k for k in keys(o.causes) if tr ∈ o.causes[k]] for tr in o.treatment]
    response_parent_names = [o.causes[k] for k in keys(o.causes) if k ∈ o.response]
    return(hcat([[intersect(tpn, rpn) for tpn in caused_by_treatment] for rpn in response_parent_names]...))
end

"""
    mediatornames(o::CausalTable, x::Symbol, y::Symbol)

Outputs the names of the mediators of the causal relationship between `x` and `y` from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the mediator names.
- `x::Symbol`, `y::Symbol`: The two variables whose mediators should be selected.

# Returns
A Vector of Symbols containing the names of the mediators between x and y.
"""
mediatornames(o::CausalTable, x::Symbol, y::Symbol) = intersect([k for k in keys(o.causes) if x ∈ o.causes[k]], o.causes[y])

"""
    mediators(o::CausalTable; collapse_parents = true)

Selects the mediators of each treatment-response pair from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the mediator variables of each treatment-response pair.
- `collape_parents::Bool`: Optional parameter, whether to collapse the output to a single `CausalTable` object if there is either only one treatment-response pair or all pair share the same set of mediators. Defaults to `true`.

# Returns
A new `CausalTable` containing only the mediators (if a single response, or all responses share the same set of mediators); otherwise, a Matrix of CausalTable objects containing the mediators of each treatment-response pair, where rows represent responses and columns represent treatments.
"""
mediators(o::CausalTable; collapse_parents = true) = select_over_matrix(o, mediatornames(o); collapse_parents = collapse_parents)

"""
    mediators(o::CausalTable, x::Symbol, y::Symbol)

Selects the mediators for a specific pair of variables (x,y) from the given `CausalTable` object; that is, the variables that are caused by x and cause y.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the mediators.
- `x::Symbol`, `y::Symbol`: The two variables whose mediators should be selected.

# Returns
A new `CausalTable` containing only the mediators of both x and y.
"""
mediators(o::CausalTable, x::Symbol, y::Symbol) = select(o, mediatornames(o, x, y))

"""
    mediatorsmatrix(o::CausalTable; collapse_parents = true)

Outputs the treatment-variable confounders from the given `CausalTable` object as a matrix (or matrix of matrices, if multiple treatment-response pairs are present).

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the confounders of each treatment-response pair.
- `collape_parents::Bool`: Optional parameter, whether to collapse the output to a single `Matrix` object if there is either only one treatment-response pair or all pair share the same set of confounders. Defaults to `true`.

# Returns
A matrix containing only the confounders.
"""
mediatorsmatrix(o::CausalTable; collapse_parents = true) = matrix(mediators(o; collapse_parents = collapse_parents))

### Instrumental Variables ###

"""
    instrumentnames(o::CausalTable)

Outputs the instrument names of each treatment-response pair from the given `CausalTable` object.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the instrument names of each treatment-response pair.

# Returns
A matrix of Vectors containing the instrument names of each treatment-response pair.
"""
#function instrumentnames(o::CausalTable)
    # Extract the mediator names: variables that cause response but are caused by treatment
#    treatment_parent_names = [o.causes[k] for k in keys(o.causes) if k ∈ o.treatment]
#    response_parent_names = [o.causes[k] for k in keys(o.causes) if k ∈ o.response]
#    return(hcat([[intersect(tpn, rpn) for tpn in caused_by_treatment] for rpn in response_parent_names]...))
#end

"""
    instrumentnames(o::CausalTable, x::Symbol, y::Symbol)

Outputs the names of the instruments of the causal relationship between `x` and `y` from the given `CausalTable` object; that is, variables that cause `x` but do not cause `y` through any other path except through `x`

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the mediator names.
- `x::Symbol`, `y::Symbol`: The two variables whose mediators should be selected.

# Returns
A Vector of Symbols containing the names of the mediators between x and y.
"""
#instrumentnames(o::CausalTable, x::Symbol, y::Symbol) = [c for c in o.causes[x]]

"""
    instruments(o::CausalTable; collapse_parents = true)

Selects the instruments of each treatment-response pair from the given `CausalTable` object; that is, variables that cause the treatment but do not cause the response through any other path except through treatment.

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the instrumental variables of each treatment-response pair.
- `collape_parents::Bool`: Optional parameter, whether to collapse the output to a single `CausalTable` object if there is either only one treatment-response pair or all pair share the same set of instruments. Defaults to `true`.

# Returns
A new `CausalTable` containing only the instruments (if a single response, or all responses share the same set of instruments); otherwise, a Matrix of CausalTable objects containing the instruments of each treatment-response pair, where rows represent responses and columns represent treatments.
"""
#instruments(o::CausalTable; collapse_parents = true) = select_over_matrix(o, instrumentnames(o); collapse_parents = collapse_parents)

"""
    instruments(o::CausalTable, x::Symbol, y::Symbol)

Selects the instruments for a specific pair of variables (x,y) from the given `CausalTable` object; that is, variables that cause `x` but do not cause `y` through any other path except through `x`

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the mediators.
- `x::Symbol`, `y::Symbol`: The two variables whose mediators should be selected.

# Returns
A new `CausalTable` containing only the mediators of both x and y.
"""
#instruments(o::CausalTable, x::Symbol, y::Symbol) = select(o, instrumentnames(o, x, y))

"""
    mediatorsmatrix(o::CausalTable; collapse_parents = true)

Outputs the treatment-variable confounders from the given `CausalTable` object as a matrix (or matrix of matrices, if multiple treatment-response pairs are present).

# Arguments
- `o::CausalTable`: The `CausalTable` object from which to extract the confounders of each treatment-response pair.
- `collape_parents::Bool`: Optional parameter, whether to collapse the output to a single `Matrix` object if there is either only one treatment-response pair or all pair share the same set of confounders. Defaults to `true`.

# Returns
A matrix containing only the confounders.
"""
#instrumentsmatrix(o::CausalTable; collapse_parents = true) = matrix(instruments(o; collapse_parents = collapse_parents))

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