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

    # Row table form
    rowtbl = CausalTable(foo3, :X, :Y, [:Z])
    @test rowtbl.tbl == foo3
    @test ncol(rowtbl) == 3
    @test nrow(rowtbl) == 3
    @test Tables.getcolumn(rowtbl, :Y) == Y

    # Column table form
    coltbl = CausalTable(foo2)
    @test coltbl.tbl == foo2
    @test ncol(coltbl) == 3
    @test nrow(coltbl) == 3
    @test Tables.columnnames(coltbl) == (:X, :Y, :Z)

    # Causal Inference
    @test gettreatment(rowtbl) == X
    @test getresponse(rowtbl) == Y
    @test getcontrols(rowtbl) == (Z = Z,)
    @test getsummaries(rowtbl) == NamedTuple()
    @test getgraph(rowtbl) == Graph()
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

    @test nrow(foo.tbl) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}
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

    dgp = DataGeneratingProcess(n -> erdos_renyi(n, 0.4), distseq; controls = [:L1, :L1_s]);
    foo = rand(dgp, 10)
    @test typeof(foo) == CausalTable
    @test Tables.columnnames(foo.tbl) == (:L1, :L1_s, :A, :A_s, :Y)

    bar = condensity(dgp, foo, :A_s)

    foo_sum = summarize(foo)

    @test nrow(foo.tbl) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}   
    @test Tables.getcolumn(foo, :L1_s) == Tables.getcolumn(summarize(foo), :L1_s)
end

