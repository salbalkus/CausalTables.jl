"""
    CausalTable

A mutable structure that contains a table (`tbl`), a graph (`graph`), and a named tuple of summaries (`summaries`).
"""
mutable struct CausalTable
    tbl
    treatment::Union{Symbol, Nothing}
    response::Union{Symbol, Nothing}
    graph::Graph
    summaries::NamedTuple
end

"""
    CausalTable(tbl)

An alternate constructor for `CausalTable` that takes a table (`tbl`) and initializes all other variables as blank.
"""
CausalTable(tbl) = CausalTable(tbl, nothing, nothing, Graph(), (;))

"""
    CausalTable(tbl)

An alternate constructor for `CausalTable` that takes a table, response, and treatment, and initializes the `graph` and `summaries` as blank.
"""
CausalTable(tbl, treatment::Union{Symbol, Nothing}, response::Union{Symbol, Nothing}) = CausalTable(tbl, treatment, response, Graph(), (;))

"""
    CausalTable(tbl)

An alternate constructor for `CausalTable` that takes a table, graph, and summaries, and initializes the `treatment` and `response` as blank.
"""
CausalTable(tbl, graph::Graph, summaries::NamedTuple) = CausalTable(tbl, nothing, nothing, graph, summaries)


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

# Basic causal inference getters
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








