using Test
using CausalTables
using Tables
using DataFrames
using Graphs
using Distributions
using Random

@testset "convolve function" begin
    @test convolve([Normal(0, 1), Normal(0, 1)]) == Normal(0, sqrt(2))
    @test convolve(Vector{Normal}(undef, 0)) == Binomial(0, 0.5)
    @test_throws ArgumentError convolve([Normal(0, 1), Uniform(0, 1)])
end

@testset "CausalTables" begin
    X = [1, 2, 3]
    Y = [4, 5, 6]
    Z = [1.0, 2.0, 3.0]
    foo1 = DataFrame(X = X, Y = Y, Z = Z)
    foo2 = Tables.columntable(foo1)
    foo3 = Tables.rowtable(foo1)
    
    # DataFrame form
    df = CausalTables.CausalTable(foo1, :X, :Y)    
    @test Tables.istable(df)
    @test Tables.columntable(foo1) == df.data
    @test ncol(df) == 3
    @test nrow(df) == 3
    @test getindex(df, 1, 2) == 4
    @test Tables.getcolumn(df, :X) == X    
    @test Tables.getcolumn(df, 1) == X
    @test Tables.columnindex(df, :X) == 1
    @test Tables.columntype(df, :X) == Int
    @test df.confounders == [:Z]

    # Row table form
    rowtbl = CausalTables.CausalTable(foo3, :X, :Y; confounders = [:Z])
    @test Tables.istable(rowtbl)
    @test Tables.rowtable(rowtbl.data) == foo3
    @test ncol(rowtbl) == 3
    @test nrow(rowtbl) == 3
    @test Tables.getcolumn(rowtbl, :Y) == Y
    @test Tables.getcolumn(rowtbl, 1) == X
    @test Tables.columnindex(rowtbl, :X) == 1
    @test Tables.columntype(rowtbl, :X) == Int

    # Column table form
    coltbl = CausalTables.CausalTable(foo2, :X, :Y; confounders = [:Z])
    @test Tables.istable(coltbl)
    @test coltbl.data == foo2
    @test ncol(coltbl) == 3
    @test nrow(coltbl) == 3
    @test Tables.getcolumn(rowtbl, :Y) == Y
    @test Tables.getcolumn(rowtbl, 1) == X
    @test Tables.columnindex(rowtbl, :X) == 1
    @test Tables.columntype(rowtbl, :X) == Int
    @test Tables.columnnames(coltbl) == (:X, :Y, :Z)

    # Extra causal-related functions
    coltbl.response
    more_sums = (S = Sum(:X, :G), T = Sum(:Z, :G), U = Sum(:Y, :G))
    coltbl2 = CausalTables.replace(coltbl, arrays = (G = [1 0 1; 0 1 1; 0 0 1],), summaries = more_sums)
    coltbl2  = summarize(coltbl2) 
    coltbl2.treatment

    @test coltbl2.treatment == [:X, :S]
    @test coltbl2.confounders == [:Z, :T]
    @test coltbl2.response == [:Y, :U]
    @test Tables.columnnames(CausalTables.treatment(coltbl2)) == (:X, :S)
    @test Tables.columnnames(CausalTables.response(coltbl2)) == (:Y, :U)
    @test Tables.columnnames(CausalTables.confounders(coltbl2)) == (:Z, :T)
    @test Tables.columnnames(CausalTables.treatmentparents(coltbl2)) == (:Z, :T)
    @test Tables.columnnames(CausalTables.responseparents(coltbl2)) == (:X, :Z, :S, :T)

    # Other convenience
    baz = (X = [4, 5], Y = ["foo", "bar"], Z = [0.1, 0.2])
    array = Graphs.adjacency_matrix(Graphs.path_graph(2))
    causalbaz = CausalTables.CausalTable(baz, :X, :Y, [:Z], (A = array,), (B = CausalTables.Sum(:X, :A),))
    @test causalbaz isa CausalTables.CausalTable
    @test CausalTables.replace(rowtbl; data = baz).data == baz
    @test Tables.subset(coltbl, 1:2).data == (X = X[1:2], Y = Y[1:2], Z = Z[1:2])
    @test Tables.subset(coltbl, 1:2; viewhint = false).data == (X = X[1:2], Y = Y[1:2], Z = Z[1:2])
    @test CausalTables.replace(rowtbl; treatment = :X).treatment == [:X]

    # Errors
    @test_throws ArgumentError CausalTables.CausalTable(foo1, :X, :X, [:Z])
    @test_throws ArgumentError CausalTables.CausalTable(foo1, :X, :Z, [:Z])
    @test_throws ArgumentError CausalTables.CausalTable(foo1, :Z, :X, [:Z])

end

@testset "DataGeneratingProcess using dgp macro, no graphs" begin
    dgp = CausalTables.@dgp(
        L1 ~ Beta(1,1),
        N = length(L1),
        L1_norm = L1 ./ sum(L1),
        L2 ~ Multinomial(N, L1_norm),
        A ~ (@. Normal(L1, 1)),
        regr = (@. A + 0.2 * L2),
        Y ~ Normal.(regr, 1)
    )

    scm = CausalTables.StructuralCausalModel(dgp, [:A], [:Y], [:L1, :L2])
    foo = rand(scm, 10)

    @test typeof(foo) == CausalTables.CausalTable
    @test Tables.columnnames(foo.data) == (:L1, :L2, :A, :Y)

    bar = CausalTables.condensity(scm, foo, :A)
    baz = CausalTables.conmean(scm, foo, :Y)
    @test nrow(foo.data) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}
    @test typeof(baz) <: Vector{T} where {T <: Real}
    @test baz == Tables.getcolumn(foo, :A) .+ 0.2 .* Tables.getcolumn(foo, :L2)

end

@testset "DataGeneratingProcess with graphs using dgp macro" begin
    dgp = @dgp(
        L1 ~ DiscreteUniform(1, 5),
        n = length(L1),
        ER = adjacency_matrix(erdos_renyi(n, 0.3)),
        L1_s $ Sum(:L1, :ER),
        A ~ (@. Normal(L1 + L1_s, 1)),
        A_s $ Sum(:A, :ER),
        Y ~ (@. Normal(A + A_s + 0.2 * L1 + 0.05 * L1_s, 1))
    )
    
    scm = CausalTables.StructuralCausalModel(dgp, [:A], [:Y], [:L1])
    foo = rand(scm, 10) 

    @test typeof(foo) == CausalTables.CausalTable
    @test Tables.columnnames(foo.data) == (:L1, :A, :Y)

    bar = CausalTables.condensity(scm, foo, :A_s)
    @test nrow(foo) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}

    foo_sum = CausalTables.summarize(foo)

    @test issetequal(Tables.columnnames(foo_sum), (:L1, :L1_s, :A, :A_s, :Y))
    @test foo.arrays.L1_s == Tables.getcolumn(foo_sum, :L1_s)
    @test foo.arrays.A_s == Tables.getcolumn(foo_sum, :A_s)

    
    # Test the graph subsetting capabilities of CausalTable
    indices = [1, 3, 7, 8]

    baz = Tables.subset(foo, indices)
    @test baz.data == Tables.subset(foo.data, indices)
    @test size(baz.arrays.ER) == (length(indices), length(indices))
end

@testset "DGP Exception throwing" begin

    ### Test that the DGP macro throws an error when it should ###

    # Test LHS
    @test_throws ArgumentError CausalTables._parse_name(:(a() ~ Normal(0, 1))) 
    @test_throws ArgumentError CausalTables._parse_name(:(1 ~ Normal(0, 1)))
    @test_throws ArgumentError CausalTables._parse_name(:(1a ~ Normal(0, 1)))
    @test_throws ArgumentError CausalTables._parse_name(:([]p23[p4] ~ Normal(0, 1)))

    expr = :(A; ~ Normal(0, 1))
    @test_throws ArgumentError CausalTables._parse_name(expr)

    ### Test the DGP constructor

    dgp = CausalTables.@dgp(
        L1 ~ DiscreteUniform(1, 5),
        A ~ (@. Normal(:L1, 1)),
        reg = :A + 0.2 * :L1,
        Y ~ (@. Normal(reg, 1))
    )

    @test_throws ArgumentError CausalTables.StructuralCausalModel(dgp; treatment = :L1, response = :Y, confounders = [:L1])
    @test_throws ArgumentError CausalTables.StructuralCausalModel(dgp; response = :L1, treatment = :A,  confounders = [:L1])
    @test_throws ArgumentError CausalTables.StructuralCausalModel(dgp; treatment = :X, response = :Y, confounders = [:L1])
    @test_throws ArgumentError CausalTables.StructuralCausalModel(dgp; treatment = :A, response = :X, confounders = [:L1])
    @test_throws ArgumentError CausalTables.StructuralCausalModel(dgp; treatment = :A, response = :Y, confounders = [:X])
    @test_throws ArgumentError CausalTables.StructuralCausalModel(dgp; treatment = :A, response = :Y, confounders = [:L1], arraynames = [:X])

    # Test the RHS
    # TODO: Currently errors in DGP construct are deferred until rand or condensity is called.
    # Can we catch them earlier?
    bad = CausalTables.StructuralCausalModel(CausalTables.@dgp(A ~ asdjfk, Y ~ adjsf); treatment = :A, response = :Y)    
    @test_throws ErrorException rand(bad, 10)
    
    tbl = CausalTables.CausalTable((A = [1, 2, 3], Y = [4, 5, 6]), treatment = :A, response = :Y)
    @test_throws ArgumentError CausalTables.condensity(bad, tbl, :not_in_dgp)
    @test_throws ErrorException CausalTables.condensity(bad, tbl, :A)
end
    
@testset "NetworkSummary" begin
    Random.seed!(1234)
    dgp = CausalTables.@dgp(
        A ~ Normal(0,1),
        L ~ Binomial(1, 0.5),
        G = adjacency_matrix(erdos_renyi(length(A), 0.3)),
        As $ Sum(:A, :G),
        Ao $ AllOrderStatistics(:A, :G),
        F $ Friends(:G),
        Lm $ Mean(:L, :G),
        Y ~ Normal(0,1)
    )
    scm = CausalTables.StructuralCausalModel(dgp; treatment = :A, response = :Y, confounders = [:L])
    tbl = rand(scm, 5)
    stbl = CausalTables.summarize(tbl)
    stbl.data.F
    @test stbl.data.As ==  stbl.arrays.G * stbl.data.A
    @test stbl.data.F == [3.0, 2.0, 2.0, 3.0, 2.0]
    @test Tables.columnnames(stbl) == (:A, :L, :Y, :As, :Ao1, :Ao2, :Ao3, :F, :Lm)
    @test stbl.treatment == [:A, :As, :Ao1, :Ao2, :Ao3]
    @test stbl.confounders == [:L, :Lm]
end