using CausalTables
using DataFrames
using Test
using Graphs
using Distributions
using Random

Random.seed!(1);

@testset "CausalTables" begin
    foo1 = DataFrame(L = [1, 2, 3], A = [5, 6, 7])
    foo2 =  (L = [1, 2, 3], A = [5, 6, 7])
    foo3 = ((L = 1, A = 5), (L = 2, A = 6), (L = 3, A = 7))

    @test CausalTable(foo1).tbl == foo1
    @test CausalTable(foo2).tbl == foo2
    @test CausalTable(foo3).tbl == foo3
end

@testset "DataGeneratingProcess" begin
    distseq = [
        :L1 => (; O...) -> DiscreteUniform(1, 5),
        :L1_s => NeighborSumIn(:L1),
        :A => (; O...) -> (@. Normal(O[:L1] + O[:L1_s], 1)),
        :A_s => NeighborSumIn(:A),
        :Y => (; O...) -> (@. Normal(O[:A] + O[:A_s] + 0.2 * O[:L1] + 0.05 * O[:L1_s], 1))
    ]

    dgp = DataGeneratingProcess(n -> erdos_renyi(n, 0.2), distseq);
    foo = rand(dgp, 10)
    
    @test typeof(foo) == CausalTable
    @test Tables.columnnames(foo.tbl) == (:L1, :L1_s, :A, :A_s, :Y)

    bar = condensity(dgp, foo, :A_s)

    @test nrow(foo.tbl) == length(bar)
    @test typeof(bar) == Vector{UnivariateDistribution}
end

