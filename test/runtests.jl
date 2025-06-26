using Test
using CausalTables
using Tables
using DataFrames
using Graphs
using Distributions
using Random
using LinearAlgebra
using DataAPI

within(x, ε) = abs(x) < ε

@testset "convolve function" begin
    @test CausalTables.convolve([Normal(0, 1), Normal(0, 1)]) == Normal(0, sqrt(2))
    @test CausalTables.convolve(Vector{Normal}(undef, 0)) == Binomial(0, 0.5)
    @test_throws ArgumentError CausalTables.convolve([Normal(0, 1), Uniform(0, 1)])
end

@testset "CausalTables" begin
    X = [1, 2, 3]
    Y = [4, 5, 6]
    Z = [1.0, 2.0, 3.0]
    foo1 = DataFrame(X = X, Y = Y, Z = Z)
    foo2 = Tables.columntable(foo1)
    foo3 = Tables.rowtable(foo1)

    # DataFrame form
    df = CausalTables.CausalTable(foo1, ["X", :Z], "Y", (X = ["Z"], Z = [], Y = ["X", :Z]))  
    @test Tables.istable(df)
    @test Tables.columntable(foo1) == df.data
    @test DataAPI.ncol(df) == 3
    @test DataAPI.nrow(df) == 3
    @test Tables.getindex(df, 1, 2) == 4
    @test Tables.getcolumn(df, :X) == X    
    @test Tables.getcolumn(df, 1) == X
    @test Tables.columnindex(df, :X) == 1
    @test Tables.columntype(df, :X) == Int
    @test df.causes == (X = [:Z], Z = [], Y = [:X, :Z])

    @test Tables.columnaccess(df)
    @test Tables.rowaccess(df)
    @test Tables.rowtable(Tables.rows(df)) == foo3
    @test Tables.schema(df) == Tables.schema(data(df))

    # Row table form
    rowtbl = CausalTables.CausalTable(foo3, :X, [:Y, :Z]; causes = (X = [], Y = [:X, :Z], Z = [:X]))
    @test Tables.istable(rowtbl)
    @test Tables.rowtable(rowtbl.data) == foo3
    @test DataAPI.ncol(rowtbl) == 3
    @test DataAPI.nrow(rowtbl) == 3
    @test Tables.getcolumn(rowtbl, :Y) == Y
    @test Tables.getcolumn(rowtbl, 1) == X
    @test Tables.columnindex(rowtbl, :X) == 1
    @test Tables.columntype(rowtbl, :X) == Int

    # Column table form
    coltbl = CausalTables.CausalTable(foo2, :X, :Y)
    @test Tables.istable(coltbl)
    @test coltbl.data == foo2
    @test DataAPI.ncol(coltbl) == 3
    @test DataAPI.nrow(coltbl) == 3
    @test Tables.getcolumn(rowtbl, :Y) == Y
    @test Tables.getcolumn(rowtbl, 1) == X
    @test Tables.columnindex(rowtbl, :X) == 1
    @test Tables.columntype(rowtbl, :X) == Int
    @test Tables.columnnames(coltbl) == (:X, :Y, :Z)

    # Extra causal-related functions
    more_sums = (S = Sum(:X, :G), T = Sum(:Z, :G), U = Sum(:Y, :G))
    coltbl2 = replace(coltbl, arrays = (G = [1 0 1; 0 1 1; 0 0 1],), summaries = more_sums)
    coltbl2 = summarize(coltbl2) 
    @test coltbl2.treatment == [:X, :S]
    coltbl2.causes

    @test coltbl2.causes == (X = [:Z, :T], Y = [:Z, :X, :T, :S], S = [:Z, :T], U = [:Z, :X, :T, :S])
    @test coltbl2.response == [:Y, :U]
    
    @test Tables.columnnames(CausalTables.treatment(coltbl2)) == (:X, :S)
    @test Tables.columnnames(CausalTables.response(coltbl2)) == (:Y, :U)
    @test Tables.columnnames(CausalTables.treatmentparents(coltbl2)) == (:Z, :T)
    @test Tables.columnnames(CausalTables.parents(coltbl2, :X)) == (:Z, :T)
    @test Tables.columnnames(CausalTables.parents(coltbl2, :Z)) == ()
    @test Tables.columnnames(CausalTables.responseparents(coltbl2)) == (:Z, :X, :T, :S)
    @test Tables.columnnames(CausalTables.parents(coltbl2, :Y)) == Tables.columnnames(CausalTables.responseparents(coltbl2))
    
    # Other convenience
    baz = (X = [4, 5], Y = ["foo", "bar"], Z = [0.1, 0.2])
    array = Graphs.adjacency_matrix(Graphs.path_graph(2))
    causalbaz = CausalTables.CausalTable(baz, :X, :Y, (X = [:Z], Y = [:X, :Z]), (A = array,), (B = CausalTables.Sum(:X, :A),))
    @test causalbaz isa CausalTables.CausalTable
    @test replace(rowtbl; data = baz).data == baz
    Tables.subset(coltbl, 1:2)
    @test Tables.subset(coltbl, 1:2).data == (X = X[1:2], Y = Y[1:2], Z = Z[1:2])
    @test Tables.subset(coltbl, 1:2; viewhint = false).data == (X = X[1:2], Y = Y[1:2], Z = Z[1:2])
    @test replace(rowtbl; treatment = :X).treatment == [:X]
    @test vec(CausalTables.treatmentmatrix(coltbl)) == X
    @test vec(CausalTables.responsematrix(coltbl)) == Y

    # Errors
    @test_throws ArgumentError CausalTables.CausalTable(foo1, :X, :X)
end

@testset "Confounders, mediators, and instrumental variables" begin
    tbl = (
        L1 = [1, 2, 3],
        L2 = [4, 5, 6],
        I1 = [0.3, 0.9, 0.7],
        I2 = [0.1, 0.2, 0.4],
        A1 = [true, false, false],
        A2 = [false, true, false],
        M1 = [true, true, false],
        M2 = [false, false, true],
        Y1 = [1.1, 2.5, 1.7],
        Y2 = [2.3, 3.4, 2.2],
        )
    causes = (A1 = [:L1, :I1], A2 = [:L1, :L2, :I1, :I2], 
         M1 = [:A1], M2 = [:A1, :A2], 
         Y1 = [:M1, :A1, :L1], 
         Y2 = [:M1, :M2, :A1, :A2, :L1, :L2])
    ctbl = CausalTable(tbl, [:A1, :A2], [:Y1, :Y2], causes)

    # Test confounders
    topcorner = hcat(tbl[:L1])
    bottomcorner = hcat(tbl[:L1], tbl[:L2])
    @test confoundersmatrix(ctbl) == Dict(:A1 => Dict(:Y1 => topcorner, :Y2 => topcorner),
                                          :A2 => Dict(:Y1 => topcorner, :Y2 => bottomcorner))
    @test confounders(ctbl, :M1, :A1).data == (;)

    # Test mediators
    @test mediatorsmatrix(ctbl) == Dict(:A1 => Dict(:Y2 => hcat(tbl[:M1], tbl[:M2]), :Y1 => hcat(tbl[:M1])),
                                        :A2 => Dict(:Y2 => hcat(tbl[:M2]), :Y1 => [;]))
    @test mediators(ctbl, :L1, :Y1).data  == (A1 = tbl[:A1],)

    # Test instrumental variables
    instrument_truth = Dict(:A1 => Dict(:Y2 => hcat(tbl[:I1]), :Y1 => hcat(tbl[:I1], tbl[:M2])),
                            :A2 => Dict(:Y2 => hcat(tbl[:I1], tbl[:I2]), :Y1 => hcat(tbl[:L2], tbl[:I1], tbl[:I2], tbl[:M2])))
    @test instrumentsmatrix(ctbl) == instrument_truth
    @test instruments(ctbl, :M1, :Y1).data == (;)

    tbl = (
        L1 = [1, 2, 3],
        I1 = [0.3, 0.9, 0.7],
        A1 = [true, false, false],
        M1 = [true, true, false],
        Y1 = [1.1, 2.5, 1.7],
        )
    causes = (A1 = [:L1, :I1], M1 = [:A1], Y1 = [:M1, :A1, :L1]) 
    ctbl = CausalTable(tbl, :A1, :Y1, causes)

    @test vec(confoundersmatrix(ctbl)) == [1,2,3]
    @test vec(mediatorsmatrix(ctbl)) == [true, true, false]
    @test vec(instrumentsmatrix(ctbl)) == [0.3, 0.9, 0.7]
end

@testset "DataGeneratingProcess utilities with no dgp macro" begin
    dgp1 = DataGeneratingProcess([O -> Normal(0,1)])
    dgp2 = DataGeneratingProcess(
        [Symbol("X$(i)") for i in 2:10],
        [O -> Normal.(O[Symbol("X$(i-1)")], 0.5) for i in 2:10]
    )
    dgp = merge(dgp1, dgp2)
    
    scm = CausalTables.StructuralCausalModel(dgp, :X1, :X2)
    ct = rand(scm, 10)

    @test typeof(ct) == CausalTables.CausalTable
    @test Tables.columnnames(ct) == Tuple(Symbol("X$(i)") for i in 1:10)
    @test_throws ArgumentError merge(DataGeneratingProcess([O -> Normal(0,1)]; varsymb = "Y"), 
                                          DataGeneratingProcess([O -> Normal(0,1)]; varsymb = "Y"))
    
    dgp3 = @dgp(
        Y ~ Normal.(reduce(+, values(O)), 1)
    )        
    scm = CausalTables.StructuralCausalModel(merge(dgp, dgp3), :X1, :Y)
    rand(scm, 10)
    @test typeof(ct) == CausalTables.CausalTable
    @test Tables.columnnames(ct) == Tuple(Symbol("X$(i)") for i in 1:10)
end

@testset "DataGeneratingProcess using dgp macro, no graphs" begin
    dgp = CausalTables.@dgp(
        L1 ~ Beta(1,1),
        N = length(L1),
        L1_norm = L1 ./ sum(L1),
        L2 ~ Multinomial(N, L1_norm),
        A ~ (@. Normal(L1, 1)),
        regr = A .+ 0.2 .* vec(sum(L2, dims=2)),
        Y ~ Normal.(regr, 1)
    )

    scm = CausalTables.StructuralCausalModel(dgp, :A, :Y)
    foo = rand(scm, 5)

    @test typeof(foo) == CausalTables.CausalTable
    @test Tables.columnnames(foo.data) == (:L1, :L2_1, :L2_2, :L2_3, :L2_4, :L2_5, :A, :Y)

    bar = CausalTables.condensity(scm, foo, :A)
    baz = CausalTables.propensity(scm, foo, :L1)
    qux = CausalTables.conmean(scm, foo, :Y)
    quux = CausalTables.convar(scm, foo, :Y)

    @test DataAPI.nrow(foo.data) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}
    @test typeof(baz) <: Vector{T} where {T <: Real}
    @test typeof(qux) <: Vector{T} where {T <: Real}
    @test typeof(quux) <: Vector{T} where {T <: Real}
    
    @test all(baz .== 1.0)
    @test qux == Tables.getcolumn(foo, :A) .+ 1.0
    @test all(quux .== 1)

    @test CausalTables.adjacency_matrix(foo) == LinearAlgebra.I
    @test CausalTables.dependency_matrix(foo) == LinearAlgebra.I

    # Test the update_arrays function
    foo_update = CausalTables.update_arrays(scm, foo)
    @test foo_update.arrays == foo.arrays
    foo2 = intervene(foo, additive_mtp(1.0))
    foo2_update = CausalTables.update_arrays(scm, foo2)
    @test all(foo2_update.arrays.regr .≈ (foo2.arrays.regr .+ 1.0))
end

@testset "DataGeneratingProcess with graphs using dgp macro" begin
    dgp = @dgp(
        L1 ~ DiscreteUniform(1, 5),
        L2 ~ DiscreteUniform(1, 5),
        n = length(L1),
        ER ≈ Graphs.adjacency_matrix(erdos_renyi(n, 0.3)),
        L1_s $ Sum(:L1, :ER),
        L2_s $ Sum(:L2, :ER),
        A ~ (@. Normal(L1 + L2 + L1_s + L2_s, 1)),
        A_s $ Sum(:A, :ER),
        μ = (@. A + A_s + 0.2 * L1 + 0.05 * L1_s),
        Y ~ (@. Normal(μ, 1))
    )
    
    scm = CausalTables.StructuralCausalModel(dgp, :A, [:Y]; causes = (A = [:L1, :L2, :L1_s], Y = [:A, :L1, :L2, :L1_s]))
    foo = rand(scm, 10) 

    @test typeof(foo) == CausalTables.CausalTable
    @test Tables.columnnames(foo.data) == (:L1, :L2, :A, :Y)

    bar = CausalTables.condensity(scm, foo, :A_s)
    @test DataAPI.nrow(foo) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}

    foo_sum = CausalTables.summarize(foo)

    @test issetequal(Tables.columnnames(foo_sum), (:L1, :L1_s, :L2, :L2_s, :A, :A_s, :Y))
    @test foo.arrays.L1_s == Tables.getcolumn(foo_sum, :L1_s)
    @test foo.arrays.L2_s == Tables.getcolumn(foo_sum, :L2_s)
    @test foo.arrays.A_s == Tables.getcolumn(foo_sum, :A_s)
    
    @test Tables.columnnames(confounders(foo)) == (:L1, :L2)
    summarize(foo).causes

    @test Tables.columnnames(confounders(summarize(foo))) == (:L1, :L2, :L1_s, :L2_s)
    @test Tables.columnnames(confounders(summarize(foo, add_summaries_as_causes = true))) == (:L1, :L2, :L1_s, :L2_s)
    @test Tables.columnnames(summarize(confounders(foo))) == (:L1, :L2, :L1_s, :L2_s)

    @test Tables.columnnames(treatment(foo)) == (:A,)
    @test Tables.columnnames(treatment(summarize(foo))) == (:A, :A_s)

    # Test the graph subsetting capabilities of CausalTable
    indices = [1, 3, 7, 8]

    baz = Tables.subset(foo, indices)
    @test baz.data == Tables.subset(foo.data, indices)
    @test size(baz.arrays.ER) == (length(indices), length(indices))

    # Test the update_arrays function
    foo_update = CausalTables.update_arrays(scm, foo)
    foo_update.arrays
    foo.arrays

    @test foo_update.arrays == foo.arrays
    foo2 = intervene(foo, additive_mtp(1.0))
    foo2_update = CausalTables.update_arrays(scm, foo2)
    @test all(foo2_update.arrays.μ .≈ (foo2.arrays.μ .+ 1 .+ sum(foo2.arrays.ER, dims=2)))
end

@testset "DGP Exception throwing" begin

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

    @test_throws ArgumentError CausalTables.StructuralCausalModel(dgp; response = :A, treatment = :A)
    @test_throws ArgumentError CausalTables.StructuralCausalModel(dgp; treatment = :X, response = :Y)
    @test_throws ArgumentError CausalTables.StructuralCausalModel(dgp; treatment = :A, response = :X)
    @test_throws ArgumentError CausalTables.StructuralCausalModel(dgp; treatment = :A, response = :Y, causes = (Y = [:A, :X],))
    @test_throws ArgumentError CausalTables.StructuralCausalModel(dgp; treatment = :A, response = :Y, arraynames = [:X])

    # Test the RHS
    # TODO: Currently errors in DGP construct are deferred until rand or condensity is called.
    # Can we catch them earlier?
    bad = CausalTables.StructuralCausalModel(CausalTables.@dgp(A ~ asdf, Y ~ dfgh); treatment = :A, response = :Y)    
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
        LoH $ AllOrderStatistics(:L, :H),
        F $ Friends(:G),
        Lm $ Mean(:L, :H),
        Y ~ Normal(0,1)
    )
    scm = CausalTables.StructuralCausalModel(dgp; treatment = :A, response = :Y)
    tbl = rand(scm, 5)

    stbl = CausalTables.summarize(tbl)

    @test stbl.data.As ==  stbl.arrays.G * stbl.data.A
    @test stbl.data.F == [2.0, 2.0, 3.0, 3.0, 2.0]
    @test Tables.columnnames(stbl) == (:A, :L, :Y, :As, :Lo1, :Lo2, :LoH1, :LoH2, :LoH3, :F, :Lm)
    @test stbl.treatment == [:A, :As]
    @test stbl.causes == (A = [:L, :Lo1, :Lo2, :LoH1, :LoH2, :LoH3, :Lm], Y = [:L, :A, :Lo1, :Lo2, :LoH1, :LoH2, :LoH3, :Lm, :As], As = [:L, :Lo1, :Lo2, :LoH1, :LoH2, :LoH3, :Lm])
    
    sub = Tables.subset(stbl, 1:3)
    @test DataAPI.nrow(sub) == 3
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
        L ~ Beta(1, 1),
        A ~ @.(Bernoulli(0.6 * L + 0.2)),
        Y ~ @.(Normal(A + 2 * L + 1))
    )

    scm = CausalTables.StructuralCausalModel(dgp, [:A], :Y)

    # Check intervention functions
    ct = rand(scm, 10000)
    cta = intervene(ct, treat_all)
    @test Tables.columnnames(cta) == Tables.columnnames(ct)
    @test all(cta.data.A .== 1.0)

    ctn = intervene(ct, treat_none)
    @test all(ctn.data.A .== 0.0)

    ε = 0.05

    # ATE
    est_ate = ate(scm)
    @test within(est_ate.μ - 1, ε)
    @test within(est_ate.eff_bound - 4.6, ε)

    # ATT
    est_att = att(scm)
    @test within(est_att.μ - 1, ε)
    @test within(est_att.eff_bound - 6.23, ε)
    # ATU
    est_atu = atu(scm)
    @test within(est_atu.μ - 1, ε)
    @test within(est_atu.eff_bound - 6.23, ε)

    # Test continuous random variables
    dgp2 = CausalTables.@dgp(
        L ~ Beta(2, 4),
        A ~ @.(Normal(L)),
        Y ~ @.(Normal(A + 2 * L + 1))
    )

    scm2 = CausalTables.StructuralCausalModel(dgp2, [:A], :Y)
    
    # Check intervention functions
    ct2 = rand(scm2, 100)
    ct_add = intervene(ct2, additive_mtp(1.0))
    @test all(ct_add.data.A .== ct2.data.A .+ 1.0)
    ct_mul = intervene(ct2, multiplicative_mtp(2.0))
    @test all(ct_mul.data.A .== ct2.data.A .* 2.0)
    
    # Modified Treatment Policy / Average Policy Effect
    est_ape_a = ape(scm2, additive_mtp(1.0))
    @test within(est_ape_a.μ - 1, ε)
    
    est_ape_m = ape(scm2, multiplicative_mtp(1.0))
    @test within(est_ape_m.μ, ε)

    mean_a = cfmean(scm2, additive_mtp(1.0))
    @test within(mean_a.μ - 3, ε)

    # Test summarized continuous random variables
    dgp3 = CausalTables.@dgp(
        L ~ Beta(1, 1),
        G = Graphs.adjacency_matrix(erdos_renyi(length(L), 3 / length(L))),
        A ~ @.(Normal(L)),
        As $ Sum(:A, :G),
        Y ~ @.(Normal(A + As + 2 * L + 1))
    )
    scm3 = CausalTables.StructuralCausalModel(dgp3, [:A, :As], :Y)
    
    # Check intervention functions
    ct3 = rand(scm3, 100)
    ct3_add = intervene(ct3, additive_mtp(1.0))
    ct3_add
    @test all(ct3_add.data.A .== ct3.data.A .+ 1.0)
    ct3_mul = intervene(ct3, multiplicative_mtp(2.0))
    @test all(ct3_mul.data.A .== ct3.data.A .* 2.0)
    
    # Modified Treatment Policy / Average Policy Effect
    est_ape_a = ape(scm3, additive_mtp(1.0))
    @test within(est_ape_a.μ - 4, ε)
    
    est_ape_m = ape(scm3, multiplicative_mtp(1.0))
    @test within(est_ape_m.μ, ε)

    mean_a = cfmean(scm3, additive_mtp(1.0))
    @test within(mean_a.μ - 8, ε)
end

@testset "Odd edge cases" begin
    d = 5
    many_distributions = DataGeneratingProcess(
        [O -> Bernoulli(1.0) for _ in 1:d]
    )
    output_distribution = @dgp(
        A ~ Normal.(reduce(+, values(O)), 0),
        Y ~ Normal.(reduce(+, values(O)), 0)
    )

    scm = StructuralCausalModel(merge(many_distributions, output_distribution), :A, :Y)

    ct = rand(scm, 10)

    @test all(ct.data.A .== d)
    @test all(ct.data.Y .== d*2)

    @test all(condensity(scm, ct, :A) .== Normal(d, 0))
    @test all(condensity(scm, ct, :Y) .== Normal(d*2, 0))
end

# Other quality checkers
include("Aqua.jl")


