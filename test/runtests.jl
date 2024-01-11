using Test
using CausalTables
using DataFrames
using Graphs
using Distributions
using Random

Random.seed!(1);

@testset "CausalTables" begin
    X = [1, 2, 3]
    Y = ["a", "b", "c"]
    Z = [1.0, 2.0, 3.0]
    foo1 = DataFrame(X = X, Y = Y, Z = Z)
    foo2 =  (X = X, Y = Y, Z = Z)
    foo3 = Tables.rowtable(((X = 1, Y = "a", Z = 1.0), (X = 2, Y = "b", Z = 2.0), (X = 3, Y = "c", Z = 3.0)))
    # DataFrame form
    df = CausalTable(foo1)
    @test df.tbl == foo1
    @test ncol(df) == 3
    @test nrow(df) == 3
    @test getindex(df, 1, 2) == "a"
    @test Tables.getcolumn(df, :X) == X
    @test gettable(df) == foo1

    # Row table form
    rowtbl = CausalTable(foo3, :X, :Y, [:Z])
    @test rowtbl.tbl == foo3
    @test ncol(rowtbl) == 3
    @test nrow(rowtbl) == 3
    @test Tables.getcolumn(rowtbl, :Y) == Y
    @test gettable(rowtbl) == foo3


    # Column table form
    coltbl = CausalTable(foo2, :X, :Y, [:Z])
    @test coltbl.tbl == foo2
    @test ncol(coltbl) == 3
    @test nrow(coltbl) == 3
    @test Tables.columnnames(coltbl) == (:X, :Y, :Z)
    @test gettable(coltbl) == foo2


    # Causal Inference
    @test gettreatmentsymbol(rowtbl) == :X
    @test getresponsesymbol(rowtbl) == :Y
    @test getcontrolssymbols(rowtbl) == [:Z]
    @test gettreatment(rowtbl) == X
    @test getresponse(rowtbl) == Y
    @test getcontrols(rowtbl).tbl == (Z = Z,)
    @test getsummaries(rowtbl) == NamedTuple()
    @test getgraph(rowtbl) == Graph()

    # Other convenience
    baz = (F = [4, 5], G = ["foo", "bar"])
    @test gettable(CausalTables.replace(rowtbl; tbl = baz)) == baz
    @test gettable(Tables.subset(coltbl, 1:2)) == (X = X[1:2], Y = Y[1:2], Z = Z[1:2])
end

@testset "DataGeneratingProcess, no graphs" begin
    distseq = [
        :L1 => (; O...) -> DiscreteUniform(1, 5),
        :A => (; O...) -> (@. Normal(O[:L1], 1)),
        :Y => (; O...) -> (@. Normal(O[:A] + 0.2 * O[:L1], 1))
    ]

    dgp = DataGeneratingProcess(distseq);
    foo = rand(dgp, 10)
    
    @test typeof(foo) == CausalTable
    @test Tables.columnnames(foo.tbl) == (:L1, :A, :Y)

    bar = condensity(dgp, foo, :A)
    baz = conmean(dgp, foo, :Y)
    @test nrow(foo.tbl) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}
    @test typeof(baz) <: Vector{T} where {T <: Real}
    @test baz == Tables.getcolumn(foo, :A) .+ 0.2 .* Tables.getcolumn(foo, :L1)
end

@testset "DataGeneratingProcess with graphs" begin
    # TODO: Make it easier to define this type of vector
    distseq = Vector{Pair{Symbol, CausalTables.ValidDGPTypes}}([
        :L1 => (; O...) -> DiscreteUniform(1, 5),
        :L1_s => NeighborSum(:L1),
        :A => (; O...) -> (@. Normal(O[:L1] + O[:L1_s], 1)),
        :A_s => NeighborSum(:A),
        :Y => (; O...) -> (@. Normal(O[:A] + O[:A_s] + 0.2 * O[:L1] + 0.05 * O[:L1_s], 1))
    ])

    dgp = DataGeneratingProcess(n -> erdos_renyi(n, 3/n), distseq; controls = [:L1, :L1_s]);
    foo = rand(dgp, 100)
    @test typeof(foo) == CausalTable
    @test Tables.columnnames(foo.tbl) == (:L1, :L1_s, :A, :A_s, :Y)

    bar = condensity(dgp, foo, :A_s)
    foo_sum = summarize(foo)
    @test nrow(foo.tbl) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}   
    @test Tables.getcolumn(foo, :L1_s) == Tables.getcolumn(summarize(foo), :L1_s)
end

