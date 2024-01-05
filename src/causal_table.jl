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
    CausalTable(tbl, treatment, response, controls, graph, summaries) = !isnothing(controls) && (treatment ∈ controls || response ∈ controls) ? error("Treatment and/or response cannot be the same as controls.") : new(tbl, treatment, response, controls, graph, summaries)
end

CausalTable(tbl) = CausalTable(tbl, nothing, nothing, nothing, Graph(), (;))
CausalTable(tbl, graph::Graph, summaries::NamedTuple) = CausalTable(tbl, nothing, nothing, nothing, graph, summaries)
# if controls not provided, assume all columns other than treatment and response are controls
CausalTable(tbl, treatment::Union{Symbol, Nothing}, response::Union{Symbol, Nothing}) = CausalTable(tbl, treatment, response, setdiff(Tables.columnnames(Tables.columns(tbl)), [treatment, response]), Graph(), (;))
CausalTable(tbl, treatment::Union{Symbol, Nothing}, response::Union{Symbol, Nothing}, controls::Vector{Symbol}) = CausalTable(tbl, treatment, response, controls, Graph(), (;))

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

# Additional convenience methods
Base.getindex(x::CausalTable, i, j) = Base.getindex(x.tbl, i, j)
DataAPI.ncol(x::CausalTable) = DataAPI.ncol(x.tbl)
DataAPI.nrow(x::CausalTable) = DataAPI.nrow(x.tbl)

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
gettreatment(x::CausalTable) = Tables.getcolumn(x, x.treatment)


"""
    getcontrols(x::CausalTable)

Selects the control variables from the given `CausalTable` object `x`.

# Arguments
- `x::CausalTable`: The `CausalTable` object from which to select the control variables.

# Returns
- `Table`: A new `CausalTable` object containing only the control variables.

"""
getcontrols(x::CausalTable) = TableOperations.select(x, x.controls)

"""
    getsummaries(x::CausalTable)

# Network causal inference getters
"""
    getsummaries(x::CausalTable)

Returns the tables stored in the CausalTable `x`.

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

function settreatment!(x::CausalTable, treatment::Symbol)
    x.treatment = treatment
end

function setresponse!(x::CausalTable, response::Symbol)
    x.response = response
end

function setcontrols!(x::CausalTable, controls::Vector{Symbol})
    x.controls = controls
end

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










