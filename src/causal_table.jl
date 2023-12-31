"""
    CausalTable

A mutable structure that contains a table (`tbl`), a graph (`graph`), and a named tuple of summaries (`summaries`).
"""
mutable struct CausalTable
    tbl
    graph::Graph
    summaries::NamedTuple
end

"""
    CausalTable(tbl)

An alternate constructor for `CausalTable` that takes a table (`tbl`) and initializes the `graph` and `summaries`.
"""
CausalTable(tbl) = CausalTable(tbl, Graph(), (;))

# declare that CausalTable is a table
Tables.istable(::Type{CausalTable}) = true

# column interface
Tables.columnaccess(::Type{<:CausalTable}) = true
Tables.columns(x::CausalTable) = Tables.columns(x.tbl)

# required Tables.AbstractColumns object methods
Tables.getcolumn(x::CausalTable, ::Type{T}, col::Int, nm::Symbol) where {T} = Tables.getcolumn(x.tbl, col)
Tables.getcolumn(x::CausalTable, nm::Symbol) = Tables.getcolumn(x.tbl, nm)
Tables.columnnames(x::CausalTable) = Tables.columnnames(x.tbl)

# fixing StackOverflow error with column indexing methods via overloading
Tables.columnindex(x::CausalTable, nm::Symbol) = Tables.columnindex(x.tbl, nm)
Tables.columntype(x::CausalTable, nm::Symbol) = Tables.columntype(x.tbl, nm)

# Additional convenience methods
Base.getindex(x::CausalTable, i) = Base.getindex(x.tbl, i)
DataAPI.ncol(x::CausalTable) = DataAPI.ncol(x.tbl)
DataAPI.nrow(x::CausalTable) = DataAPI.nrow(x.tbl)






