using Test
using CausalTables
using Tables
using DataFrames
using Graphs
using Distributions
using Random
using LinearAlgebra

within(x, ε) = abs(x) < ε

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
    @test Tables.getindex(df, 1, 2) == 4
    @test Tables.getcolumn(df, :X) == X    
    @test Tables.getcolumn(df, 1) == X
    @test Tables.columnindex(df, :X) == 1
    @test Tables.columntype(df, :X) == Int
    @test df.confounders == [:Z]

    @test Tables.columnaccess(df)
    @test Tables.rowaccess(df)
    @test Tables.rowtable(Tables.rows(df)) == foo3
    @test Tables.schema(df) == Tables.schema(data(df))

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

    @test coltbl2.treatment == [:X, :S]
    @test coltbl2.confounders == [:Z, :T]
    @test coltbl2.response == [:Y, :U]
    @test Tables.columnnames(CausalTables.treatment(coltbl2)) == (:X, :S)
    @test Tables.columnnames(CausalTables.response(coltbl2)) == (:Y, :U)
    @test Tables.columnnames(CausalTables.confounders(coltbl2)) == (:Z, :T)
    @test Tables.columnnames(CausalTables.treatmentparents(coltbl2)) == (:Z, :T)
    @test Tables.columnnames(CausalTables.parents(coltbl2, :X)) == (:Z, :T)
    @test Tables.columnnames(CausalTables.parents(coltbl2, :Z)) == ()
    @test Tables.columnnames(CausalTables.responseparents(coltbl2)) == (:X, :Z, :S, :T)
    @test Tables.columnnames(CausalTables.parents(coltbl2, :Y)) == Tables.columnnames(CausalTables.responseparents(coltbl2))

    # Other convenience
    baz = (X = [4, 5], Y = ["foo", "bar"], Z = [0.1, 0.2])
    array = Graphs.adjacency_matrix(Graphs.path_graph(2))
    causalbaz = CausalTables.CausalTable(baz, :X, :Y, [:Z], (A = array,), (B = CausalTables.Sum(:X, :A),))
    @test causalbaz isa CausalTables.CausalTable
    @test CausalTables.replace(rowtbl; data = baz).data == baz
    @test Tables.subset(coltbl, 1:2).data == (X = X[1:2], Y = Y[1:2], Z = Z[1:2])
    @test Tables.subset(coltbl, 1:2; viewhint = false).data == (X = X[1:2], Y = Y[1:2], Z = Z[1:2])
    @test CausalTables.replace(rowtbl; treatment = :X).treatment == [:X]
    @test vec(CausalTables.treatmentmatrix(coltbl)) == X
    @test vec(CausalTables.responsematrix(coltbl)) == Y
    @test vec(CausalTables.confoundersmatrix(coltbl)) == Z

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

    scm = CausalTables.StructuralCausalModel(dgp, :A, :Y, [:L1, :L2])
    foo = rand(scm, 10)

    @test typeof(foo) == CausalTables.CausalTable
    @test Tables.columnnames(foo.data) == (:L1, :L2, :A, :Y)

    bar = CausalTables.condensity(scm, foo, :A)
    baz = CausalTables.propensity(scm, foo, :L1)
    qux = CausalTables.conmean(scm, foo, :Y)
    quux = CausalTables.convar(scm, foo, :Y)

    @test nrow(foo.data) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}
    @test typeof(baz) <: Vector{T} where {T <: Real}
    @test typeof(qux) <: Vector{T} where {T <: Real}
    @test typeof(quux) <: Vector{T} where {T <: Real}
    
    @test all(baz .== 1.0)
    @test qux == Tables.getcolumn(foo, :A) .+ 0.2 .* Tables.getcolumn(foo, :L2)
    @test all(quux .== 1)

    @test CausalTables.adjacency_matrix(foo) == LinearAlgebra.I
    @test CausalTables.dependency_matrix(foo) == LinearAlgebra.I
end

@testset "DataGeneratingProcess with graphs using dgp macro" begin
    dgp = @dgp(
        L1 ~ DiscreteUniform(1, 5),
        n = length(L1),
        ER = Graphs.adjacency_matrix(erdos_renyi(n, 0.3)),
        L1_s $ Sum(:L1, :ER),
        A ~ (@. Normal(L1 + L1_s, 1)),
        A_s $ Sum(:A, :ER),
        Y ~ (@. Normal(A + A_s + 0.2 * L1 + 0.05 * L1_s, 1))
    )
    
    scm = CausalTables.StructuralCausalModel(dgp, :A, [:Y], [:L1])
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
        L ~ Binomial(6, 0.5),
        G = Graphs.adjacency_matrix(erdos_renyi(length(L), 0.5)),
        H = Graphs.adjacency_matrix(erdos_renyi(length(L), 0.5)),
        As $ Sum(:A, :G),
        Lo $ KOrderStatistics(:L, :G, 2),
        F $ Friends(:G),
        Lm $ Mean(:L, :H),
        Y ~ Normal(0,1)
    )
    scm = CausalTables.StructuralCausalModel(dgp; treatment = :A, response = :Y, confounders = [:L])
    tbl = rand(scm, 5)

    stbl = CausalTables.summarize(tbl)

    @test stbl.data.As ==  stbl.arrays.G * stbl.data.A
    @test stbl.data.F == [2.0, 2.0, 3.0, 3.0, 2.0]
    @test Tables.columnnames(stbl) == (:A, :L, :Y, :As, :Lo1, :Lo2, :F, :Lm)
    @test stbl.treatment == [:A, :As]
    @test stbl.confounders == [:L, :Lo1, :Lo2, :Lm]
    
    sub = Tables.subset(stbl, 1:3)
    @test nrow(sub) == 3
    @test size(sub.arrays.G) == (3, 3)

    adj = CausalTables.adjacency_matrix(tbl)
    @test sum(adj) == 18
    @test all(map(x -> x ∈ [0.0, 1.0], vec(adj)))

    dep = CausalTables.dependency_matrix(tbl)
    @test sum(dep) == 25
    @test all(map(x -> x ∈ [0.0, 1.0], vec(dep)))
end

@testset "Counterfactual estimand approximation" begin
    Random.seed!(1234)

    # Test binary random variables
    dgp = CausalTables.@dgp(
        L ~ Beta(2, 4),
        A ~ @.(Bernoulli(L)),
        Y ~ @.(Normal(A + 2 * L + 1))
    )

    scm = CausalTables.StructuralCausalModel(dgp, [:A], :Y)

    # Check intervention functions
    ct = rand(scm, 100)
    cta = intervene(ct, treat_all)
    @test Tables.columnnames(cta) == Tables.columnnames(ct)
    @test all(cta.data.A .== 1.0)

    ctn = intervene(ct, treat_none)
    @test all(ctn.data.A .== 0.0)

    ε = 0.05
    # ATE
    est_ate = ate(scm)
    @test within(est_ate.μ - 1, ε)
    @test within(est_ate.eff_bound - 2, ε)

    # ATT
    est_att = att(scm)
    est_att.μ - 1
    @test within(est_att.μ - 1, ε)
    @test within(est_att.eff_bound - 2, ε)
    # ATU
    est_atu = atu(scm)
    @test within(est_atu.μ - 1, ε)
    @test within(est_atu.eff_bound - 2, ε)

    # Test continuous random variables
    dgp2 = CausalTables.@dgp(
        L ~ Beta(2, 4),
        A ~ @.(Normal(L)),
        Y ~ @.(Normal(A + 2 * L + 1))
    )

    scm2 = CausalTables.StructuralCausalModel(dgp2, [:A], :Y, [:L])
    
    # Check intervention functions
    ct2 = rand(scm2, 100)
    ct_add = intervene(ct2, additive_mtp(1.0))
    @test all(ct_add.data.A .== ct2.data.A .+ 1.0)
    ct_mul = intervene(ct2, multiplicative_mtp(2.0))
    @test all(ct_mul.data.A .== ct2.data.A .* 2.0)
    
    # Modified Treatment Policy / Average Policy Effect
    est_ape_a = ape(scm2, additive_mtp(1.0))
    @test within(est_ape_a.μ - 1, ε)
    @test within(est_ape_a.eff_bound - 2, ε)
    
    est_ape_m = ape(scm2, multiplicative_mtp(1.0))
    @test within(est_ape_m.μ, ε)
    @test within(est_ape_m.eff_bound - 2, ε)

    mean_a = cfmean(scm2, additive_mtp(1.0))
    @test within(mean_a.μ - 3, ε)
    @test within(mean_a.eff_bound - 2.3, 0.1)
end

