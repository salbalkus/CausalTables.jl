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
    length(intersect(treatment, response, confounders)) != 0 && throw(ArgumentError("Treatment, response, and confounder sets must be disjoint")) 
    
    # Return fully processed causal variable names
    return treatment, response, confounders
end

mutable struct CausalTable

    # data storage
    data::Tables.MatrixTable

    # labels
    treatment::Symbols
    response::Symbols
    confounders::Symbols

    # other
    arrays::NamedTuple
    summaries::NamedTuple

    function CausalTable(data, treatment, response, confounders, arrays, summaries)
        
        ## Process treatment and response variables into vectors, if they are not already vectors
        treatment, response, confounders = _process_causal_variable_names(treatment, response, confounders)

        # Ensure data input is a Table
        !Tables.istable(data) && throw(ArgumentError("`data` must be a Table. See https://tables.juliadata.org/ for more information."))        

        # Ensure treatment, response, and confounders are contained within the data
        names = columnnames(columns(data))
        any(t ∉ names for t in treatment)   && throw(ArgumentError("Treatment variable(s) not found in data"))
        any(r ∉ names for r in response)    && throw(ArgumentError("Response variable(s) not found in data"))
        any(c ∉ names for c in confounders) && throw(ArgumentError("Confounder variable(s) not found in data"))

        ## Construction ##

        # store the names of the input data columns
        names = Tables.columnnames(Tables.columns(data))  

        # store a matrix of data from the input table  
        data_matrix = Tables.matrix(data)
        
        # Construct a CausalTable with an underlying MatrixTable to store random vectors
        new(Tables.table(data_matrix; header = names), treatment, response, confounders, arrays, summaries)
    end
end

CausalTable(data, treatment, response; confounders = [], arrays = (;), summaries = (;)) = CausalTable(data, treatment, response, confounders, arrays, summaries)

Tables.istable(::Type{CausalTable}) = true

### Column Interface ###
# Currently only allow column access from fixed data table (not network)

Tables.columnaccess(::Type{CausalTable}) = true
Tables.columns(o::CausalTable) = Tables.columns(o.data)

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
    data_subset = Tables.subset(o.data, inds; viewhint)

    if isnothing(viewhint) || viewhint
        arrays_subset = map(x -> view(x, inds, inds), o.arrays)
    else
        arrays_subset = map(x -> getindex(x, inds, inds), o.arrays)
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
replace(o::CausalTable; kwargs...) = CausalTable([field in keys(kwargs) ?  kwargs[field] : getfield(O, field) for field in fieldnames(typeof(O))]...)


getscm(o::CausalTable) = merge(Tables.columntable(o.data), o.arrays)