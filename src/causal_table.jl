"""
    CausalTable

A mutable structure that contains a table (`tbl`), a graph (`graph`), and a named tuple of summaries (`summaries`).
"""
mutable struct CausalTable
    tbl
    treatment::Union{Symbol, Nothing}
    response::Union{Symbol, Nothing}
    controls::Union{Vector{Symbol}, Nothing}
    graph::Graph
    summaries::NamedTuple
    function CausalTable(tbl, treatment, response, controls, graph, summaries)
        cols = Tables.columns(tbl)
        if !isnothing(controls) && (treatment ∈ controls || response ∈ controls) 
            throw(ArgumentError("Treatment and/or response cannot be the same as controls.")) 
        elseif !isnothing(treatment) && !isnothing(response) && treatment == response
            throw(ArgumentError("Treatment cannot be the same as response"))
        else 
            new(tbl, treatment, response, controls, graph, summaries)
        end
    end
end

CausalTable(tbl) = CausalTable(tbl, nothing, nothing, nothing, Graph(), (;))
CausalTable(tbl, graph::Graph, summaries::NamedTuple) = CausalTable(tbl, nothing, nothing, nothing, graph, summaries)

# if controls not provided, assume all columns other than treatment and response are controls
CausalTable(tbl, treatment::Union{Symbol, Nothing}, response::Union{Symbol, Nothing}) = CausalTable(tbl, treatment, response, setdiff(Tables.columnnames(Tables.columns(tbl)), [treatment, response]), Graph(), (;))
CausalTable(tbl, treatment::Union{Symbol, Nothing}, response::Union{Symbol, Nothing}, controls::Vector{Symbol}) = CausalTable(tbl, treatment, response, controls, Graph(), (;))

# Construct using kwargs
CausalTable(tbl = nothing; treatment = nothing, response = nothing, controls = nothing, graph = Graph(), summaries = (;)) = CausalTable(tbl, treatment, response, controls, graph, summaries)

# declare that CausalTable is a table
Tables.istable(::Type{CausalTable}) = true

# column interface
Tables.columnaccess(::Type{<:CausalTable}) = true
Tables.columns(x::CausalTable) = Tables.columns(x.tbl)

# TODO: Is casting to Tables.columns too slow?
# required Tables.AbstractColumns object methods
Tables.getcolumn(x::CausalTable, ::Type{T}, col::Int, nm::Symbol) where {T} = Tables.getcolumn(Tables.columns(x.tbl), col)
Tables.getcolumn(x::CausalTable, nm::Symbol) = Tables.getcolumn(Tables.columns(x.tbl), nm)
Tables.columnnames(x::CausalTable) = Tables.columnnames(Tables.columns(x.tbl))

# fixing StackOverflow error with column indexing methods via overloading
Tables.columnindex(x::CausalTable, nm::Symbol) = Tables.columnindex(x.tbl, nm)
Tables.columntype(x::CausalTable, nm::Symbol) = Tables.columntype(x.tbl, nm)

# Additional convenience functions from existing method
DataAPI.ncol(x::CausalTable) = DataAPI.ncol(x.tbl)
DataAPI.nrow(x::CausalTable) = DataAPI.nrow(x.tbl)
Base.getindex(x::CausalTable, i::Int, j::Int) = Base.getindex(x.tbl, i, j)

# Basic causal inference getters and setters
"""
    getresponse(x::CausalTable)

Get the response variable from a CausalTable.

# Arguments
- `x::CausalTable`: The CausalTable object.

# Returns
The response variable from the CausalTable.
"""
getresponse(x::CausalTable) = Tables.getcolumn(x, x.response)
"""
    gettreatment(x::CausalTable)

Get the treatment column from a CausalTable.

# Arguments
- `x::CausalTable`: The CausalTable object.

# Returns
The treatment column of the CausalTable.
"""
gettreatment(x::CausalTable) = Tables.getcolumn(x, gettreatmentsymbol(x))

"""
    getcontrols(x::CausalTable)

Selects the control variables from the given `CausalTable` object `x`.

# Arguments
- `x::CausalTable`: The `CausalTable` object from which to select the control variables.
- `keepcausal::Bool`: Determines whether to keep the CausalTable wrapping or return a NamedTuple. Default is `true`.

# Returns
A new `CausalTable` object containing only the control variables.

"""
function getcontrols(x::CausalTable; keepcausal = true)
    L = TableOperations.select(x, x.controls...) |> Tables.columntable
    if keepcausal
        L = CausalTable(L, gettreatmentsymbol(x), 
                           getresponsesymbol(x), 
                           getcontrolssymbols(x), 
                           getgraph(x), 
                           getsummaries(x))
    end
    return L
end

gettreatmentsymbol(x::CausalTable) = x.treatment
getresponsesymbol(x::CausalTable) = x.response
getcontrolssymbols(x::CausalTable) = x.controls


# Network causal inference getters
"""
    getsummaries(x::CausalTable)

# Arguments
- `x::CausalTable`: The CausalTable object.

# Returns
An array of tables stored in the CausalTable `x`.
"""
getsummaries(x::CausalTable) = x.summaries
"""
    getgraph(x::CausalTable)

Get the graph associated with a CausalTable.

# Arguments
- `x::CausalTable`: The CausalTable object.

# Returns
- The graph associated with the CausalTable.
"""
getgraph(x::CausalTable) = x.graph


"""
    gettable(x::CausalTable)

Extracts the underlying table from a `CausalTable`.

# Arguments
- `x::CausalTable`: The `CausalTable` object.

# Returns
- The underlying table.

"""
gettable(x::CausalTable) = x.tbl


"""
    settreatment!(x::CausalTable, treatment::Symbol)

Set the treatment variable for a CausalTable.

# Arguments
- `x::CausalTable`: The CausalTable object.
- `treatment::Symbol`: The symbol representing the new treatment variable.

"""
function settreatment!(x::CausalTable, treatment::Symbol)
    x.treatment = treatment
end

"""
    setresponse!(x::CausalTable, response::Symbol)

Set the response variable for a CausalTable.

# Arguments
- `x::CausalTable`: The CausalTable object.
- `response::Symbol`: The symbol representing the new response variable.

"""
function setresponse!(x::CausalTable, response::Symbol)
    x.response = response
end

"""
    setcontrols!(x::CausalTable, controls::Vector{Symbol})

Set the control variables for a CausalTable.

# Arguments
- `x::CausalTable`: The CausalTable object.
- `controls::Vector{Symbol}`: The new control variables to be set.

"""
function setcontrols!(x::CausalTable, controls::Vector{Symbol})
    x.controls = controls
end

"""
    setcausalvars!(x::CausalTable; treatment=nothing, response=nothing, controls=nothing)

Convenience function for setting new treatment, response, and controls variables for a CausalTable at once.

Arguments:
- `x::CausalTable`: The CausalTable object.
- `treatment=nothing`: The treatment variable.
- `response=nothing`: The response variable.
- `controls=nothing`: The control variables.
"""
function setcausalvars!(x::CausalTable; treatment=nothing, response=nothing, controls=nothing)
    if !isnothing(treatment)
        settreatment!(x, treatment)
    end
    if !isnothing(response)
        setresponse!(x, response)
    end
    if !isnothing(controls)
        setcontrols!(x, controls)
    end
end

"""
    replace(x::CausalTable; tbl = nothing, treatment = nothing, response = nothing, controls = nothing, graph = nothing, summaries = nothing)

Conviently replace several components of a CausalTable with new values.

Arguments:
- `x::CausalTable`: The CausalTable object to modify.
- `tbl`: The new table to replace the existing table. If `nothing`, the current table is used.
- `treatment`: The new treatment symbol to replace the existing treatment symbol. If `nothing`, the current treatment symbol is used.
- `response`: The new response symbol to replace the existing response symbol. If `nothing`, the current response symbol is used.
- `controls`: The new control symbols to replace the existing control symbols. If `nothing`, the current control symbols are used.
- `graph`: The new graph to replace the existing graph. If `nothing`, the current graph is used.
- `summaries`: The new summaries to replace the existing summaries. If `nothing`, the current summaries are used.

Returns:
- `CausalTable`: The modified CausalTable object.

"""
function replace(x::CausalTable; tbl = nothing, treatment = nothing, response = nothing, controls = nothing, graph = nothing, summaries = nothing)
    if isnothing(tbl)
        tbl = gettable(x)
    end
    if isnothing(treatment)
        treatment = gettreatmentsymbol(x)
    end
    if isnothing(response)
        response = getresponsesymbol(x)
    end
    if isnothing(controls)
        controls = getcontrolssymbols(x)
    end
    if isnothing(graph)
        graph = getgraph(x)
    end
    if isnothing(summaries)
        summaries = getsummaries(x)
    end

    return CausalTable(tbl, treatment, response, controls, graph, summaries)
end

"""
    replacetable(x::CausalTable, tbl)

Conveniently replace the underlying table of a `CausalTable` with a new table.

# Arguments
- `x::CausalTable`: The `CausalTable` object to replace the table for.
- `tbl`: The new table to replace the existing table with.

# Returns
A new `CausalTable` object with the updated table.

"""
replacetable(x::CausalTable, tbl) = CausalTable(tbl, x.treatment, x.response, x.controls, x.graph, x.summaries)


"""
    subset(x::CausalTable, ind)

Subset a CausalTable `x` based on the given indices `ind`. Note that viewhinting is not supported; this function will return a *copy* of the CausalTable.

# Arguments
- `x`: The CausalTable to be subsetted.
- `ind`: The indices to subset the table and graph.

# Returns
A new CausalTable with the subsetted table. If the `graph` attribute is not Nothing, then this function takes the subgraph induced by the subsetted indices using `induced`

"""
function Tables.subset(x::CausalTable, ind; viewhint = nothing)
    # First check to see if user is trying to viewhint; Graphs do not support this, so we need to throw an error
    if viewhint == true
        throw(ArgumentError("CausalTables do not support viewhinting. Remove `viewhint = true` argument from your call to Tables.subset"))
    end

    graph_subset = getgraph(x) # default graph, for when graph has no edges

    # Subset the table
    tbl_subset = Tables.subset(gettable(x), ind)

    # Subset the graph
    if nv(graph_subset) > 0 # If the graph isn't empty...
        if length(unique(ind)) != length(ind) # If there are duplicate indices...
            throw(ArgumentError("Cannot subset CausalTable with non-unique indices."))
        else
            graph_subset = induced_subgraph(graph_subset, ind)[1] # Subset the graph
        end
    end

    return replace(x; tbl = tbl_subset, graph = graph_subset)
end













