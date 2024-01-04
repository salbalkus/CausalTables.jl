using CausalTables
using DataFrames
using Test
using Graphs
using Distributions
using Random

Random.seed!(1);

@testset "CausalTables" begin
    X = [1, 2, 3]
    Y = ["a", "b", "c"]
    foo1 = DataFrame(X = X, Y = Y)
    foo2 =  (X = X, Y = Y)
    foo3 = Tables.rowtable(((X = 1, Y = "a"), (X = 2, Y = "b"), (X = 3, Y = "c")))

    # DataFrame form
    df = CausalTable(foo1)
    @test df.tbl == foo1
    @test ncol(df) == 2
    @test nrow(df) == 3
    @test getindex(df, 1, 2) == "a"
    @test Tables.getcolumn(df, :X) == X

    # Row table form
    rowtbl = CausalTable(foo3, :X, :Y)
    @test rowtbl.tbl == foo3
    @test ncol(rowtbl) == 2
    @test nrow(rowtbl) == 3
    @test Tables.getcolumn(rowtbl, :Y) == Y

    # Column table form
    coltbl = CausalTable(foo2)
    @test coltbl.tbl == foo2
    @test ncol(coltbl) == 2
    @test nrow(coltbl) == 3
    @test Tables.columnnames(coltbl) == (:X, :Y)

    # Causal Inference
    @test gettreatment(rowtbl) == X
    @test getresponse(rowtbl) == Y
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
        :L1_s => NeighborSumIn(:L1),
        :A => (; O...) -> (@. Normal(O[:L1] + O[:L1_s], 1)),
        :A_s => NeighborSumIn(:A),
        :Y => (; O...) -> (@. Normal(O[:A] + O[:A_s] + 0.2 * O[:L1] + 0.05 * O[:L1_s], 1))
    ])

    dgp = DataGeneratingProcess(n -> erdos_renyi(n, 0.4), distseq);
    foo = rand(dgp, 10)
    Tables.getcolumn(foo, :L1_s)
    @test typeof(foo) == CausalTable
    @test Tables.columnnames(foo.tbl) == (:L1, :L1_s, :A, :A_s, :Y)

    bar = condensity(dgp, foo, :A_s)

    foo_sum = summarize(foo)

    @test nrow(foo.tbl) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}   
    @test Tables.getcolumn(foo, :L1_s) == Tables.getcolumn(summarize(foo), :L1_s)
end

