# Define a new type, a Vector of Symbols

# Define Exception messages
SUMMARY_NOTE = "Note: If response is a summary over a network (contained within tbl.summaries), make sure that you call `summary(tbl::CausalTable)` on your table before calling"

"""
    CausalTable

A mutable structure that contains a table (`tbl`), a graph (`graph`), and a named tuple of summaries (`summaries`).
"""
mutable struct CausalTable
    tbl
    treatment::NamedTuple
    response::NamedTuple
    control::SymbolVector
    net::SymbolVector
    summaries::NamedTuple
    function CausalTable(tbl, treatment, response, control, net, summaries)
        if length(intersect(treatment, response, control, net)) != 0
            throw(ArgumentError("Treatment, response, controls, and networks must be different variables")) 
            new(tbl, treatment, response, control, net, summaries)
        end
    end
end

# Construct using kwargs
CausalTable(tbl = nothing; treatment = [], response = [], control = [], net = [], summaries = (;)) = CausalTable(tbl, treatment, response, control, net, summaries)

# declare that CausalTable is a table
Tables.istable(::Type{CausalTable}) = true

# column interface
Tables.columnaccess(::Type{<:CausalTable}) = true
Tables.columns(x::CausalTable) = Tables.columns(x.tbl)

# TODO: Is casting to Tables.columns too slow?
# required Tables.AbstractColumns object methods
Tables.getcolumn(x::CausalTable, nm::Symbol) = Tables.getcolumn(Tables.columns(x.tbl), nm)
Tables.getcolumn(x::CausalTable, col::Int) = Tables.getcolumn(Tables.columns(x.tbl), col)
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

Get the response from a CausalTable.

# Arguments
- `x::CausalTable`: The CausalTable object.

# Returns
The response variable(s) from the CausalTable.
"""
function getresponse(x::CausalTable; keepcausal = true)
    try
        L = TableOperations.select(x, getresponsesymbol(x)...) |> Tables.columntable
    catch
        error("One or more response variables not contained in the data. $(SUMMARY_NOTE) `getcontrols`.")
    end 
    if keepcausal
        L = CausalTable(L, gettreatmentsymbol(x), 
                        getresponsesymbol(x), 
                        getcontrolsymbol(x), 
                        getnetsymbol(x), 
                        getsummaries(x))
    end
    return L
end

"""
    gettreatment(x::CausalTable)

Get the treatment from a CausalTable.

# Arguments
- `x::CausalTable`: The CausalTable object.

# Returns
The treatment variable(s) of the CausalTable.
"""
function gettreatment(x::CausalTable; keepcausal = true)
    try
        L = TableOperations.select(x, gettreatmentsymbol(x)...) |> Tables.columntable
    catch
        error("One or more treatment variables not contained in the data. $(SUMMARY_NOTE) `gettreatment`.")
    end 
    if keepcausal
        L = CausalTable(L, gettreatmentsymbol(x), 
                        getresponsesymbol(x), 
                        getcontrolsymbol(x), 
                        getnetsymbol(x), 
                        getsummaries(x))
    end
    return L
end
"""
    getcontrol(x::CausalTable)

Selects the control variables from the given `CausalTable` object `x`.

# Arguments
- `x::CausalTable`: The `CausalTable` object from which to select the control variables.
- `keepcausal::Bool`: Determines whether to keep the CausalTable wrapping or return a NamedTuple. Default is `true`.

# Returns
A new `CausalTable` object containing only the control variables.

"""
function getcontrol(x::CausalTable; keepcausal = true)
    try
        L = TableOperations.select(x, x.controls...) |> Tables.columntable
    catch
        error("One or more control variables not contained in the data. $(SUMMARY_NOTE) `getcontrol`.")
    end 
    if keepcausal
        L = CausalTable(L, gettreatmentsymbol(x), 
                        getresponsesymbol(x), 
                        getcontrolsymbol(x), 
                        getnetsymbol(x), 
                        getsummaries(x))
    end
    return L
end

gettreatmentsymbol(x::CausalTable) = x.treatment
getresponsesymbol(x::CausalTable) = x.response
getcontrolsymbol(x::CausalTable) = x.control
getnetsymbol(x::CausalTable) = x.net


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
    getnet(x::CausalTable)

Get the network(s) associated with a CausalTable.

# Arguments
- `x::CausalTable`: The CausalTable object.

# Returns
- The network(s) associated with the CausalTable.
"""
function getnet(x::CausalTable; keepcausal = true)
    try
        L = TableOperations.select(x, getnetsymbol(x)...) |> Tables.columntable
    catch
        error("One or more control variables not contained in the data. $(SUMMARY_NOTE) `getcontrol`.")
    end 
    if keepcausal
        L = CausalTable(L, gettreatmentsymbol(x), 
                        getresponsesymbol(x), 
                        getcontrolsymbol(x), 
                        getnetsymbol(x), 
                        getsummaries(x))
    end
    return L
end


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
function settreatment!(x::CausalTable, treatment::SymbolVector)
    x.treatment = treatment
end

"""
    setresponse!(x::CausalTable, response::Symbol)

Set the response variable for a CausalTable.

# Arguments
- `x::CausalTable`: The CausalTable object.
- `response::Symbol`: The symbol representing the new response variable.

"""
function setresponse!(x::CausalTable, response::SymbolVector)
    x.response = response
end

"""
    setcontrols!(x::CausalTable, controls::Vector{Symbol})

Set the control variables for a CausalTable.

# Arguments
- `x::CausalTable`: The CausalTable object.
- `controls::Vector{Symbol}`: The new control variables to be set.

"""
function setcontrol!(x::CausalTable, control::SymbolVector)
    x.control = control
end

function setnet!(x::CausalTable, net::SymbolVector)
    x.net = net
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
function setcausalvars!(x::CausalTable; treatment=nothing, response=nothing, control=nothing, net=nothing)
    if !isnothing(treatment)
        settreatment!(x, treatment)
    end
    if !isnothing(response)
        setresponse!(x, response)
    end
    if !isnothing(control)
        setcontrols!(x, controls)
    end
    if !isnothing(net)
        setnet!(x, net)
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
function replace(x::CausalTable; tbl = nothing, treatment = nothing, response = nothing, control = nothing, graph = nothing, summaries = nothing)
    if isnothing(tbl)
        tbl = gettable(x)
    end
    if isnothing(treatment)
        treatment = gettreatmentsymbol(x)
    end
    if isnothing(response)
        response = getresponsesymbol(x)
    end
    if isnothing(control)
        control = getcontrolsymbol(x)
    end
    if isnothing(net)
        net = getnetsymbol(x)
    end
    if isnothing(summaries)
        summaries = getsummaries(x)
    end

    return CausalTable(tbl, treatment, response, control, net, summaries)
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
replacetable(x::CausalTable, tbl) = CausalTable(tbl, x.treatment, x.response, x.controls, x.net, x.summaries)


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

    net_subset = 1:n # default graph, for when graph has no edges

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













