using Test
using CausalTables
using Tables
using TableOperations
import StatsBase: mode
using DataFrames
using Graphs
using Distributions
using Random
using LinearAlgebra

@testset "convolve function" begin
    @test convolve([Normal(0, 1), Normal(0, 1)]) == Normal(0, sqrt(2))
    @test_throws ArgumentError convolve(Vector{Normal}(undef, 0))
    @test_throws ArgumentError convolve([Normal(0, 1), Uniform(0, 1)])
end

a = [:A => 1, :B => 2, :C => 3]
tag = [:data, :data, :arrays]
a[tag .== :data]

m = !=(0).([1 0 1; 0 1 0; 1 0 1])
using Distributions
a = [Normal(1, 1), Normal(2, 1), Normal(3, 1)]



@testset "CausalTables" begin
    X = [1, 2, 3]
    Y = ["a", "b", "c"]
    Z = [1.0, 2.0, 3.0]
    foo1 = DataFrame(X = X, Y = Y, Z = Z)
    foo2 =  (X = X, Y = Y, Z = Z)
    foo3 = Tables.rowtable(((X = 1, Y = "a", Z = 1.0), (X = 2, Y = "b", Z = 2.0), (X = 3, Y = "c", Z = 3.0)))
    
    # DataFrame form
    df = CausalTable(foo1)

    @test Tables.istable(df)
    @test df.tbl == foo1
    @test ncol(df) == 3
    @test nrow(df) == 3
    @test getindex(df, 1, 2) == "a"
    @test Tables.getcolumn(df, :X) == X    
    @test Tables.getcolumn(df, 1) == X
    @test Tables.columnindex(df, :X) == 1
    @test Tables.columntype(df, :X) == Int
    @test gettable(df) == foo1

    # Row table form
    rowtbl = CausalTable(foo3, :X, :Y, [:Z])
    @test Tables.istable(rowtbl)
    @test rowtbl.tbl == foo3
    @test ncol(rowtbl) == 3
    @test nrow(rowtbl) == 3
    @test Tables.getcolumn(rowtbl, :Y) == Y
    @test Tables.getcolumn(rowtbl, 1) == X
    @test Tables.columnindex(rowtbl, :X) == 1
    @test Tables.columntype(rowtbl, :X) == Int
    @test gettable(rowtbl) == foo3


    # Column table form
    coltbl = CausalTable(foo2, :X, :Y, [:Z])
    @test Tables.istable(coltbl)
    @test coltbl.tbl == foo2
    @test ncol(coltbl) == 3
    @test nrow(coltbl) == 3
    @test Tables.getcolumn(rowtbl, :Y) == Y
    @test Tables.getcolumn(rowtbl, 1) == X
    @test Tables.columnindex(rowtbl, :X) == 1
    @test Tables.columntype(rowtbl, :X) == Int
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
    baz = (X = [4, 5], Y = ["foo", "bar"], Z = [0.1, 0.2])
    
    @test CausalTable(baz, path_graph(2), (X_s = Sum(:X),)) isa CausalTable
    @test CausalTable(baz, :X, :Y) isa CausalTable

    @test gettable(CausalTables.replace(rowtbl; tbl = baz)) == baz
    @test gettable(Tables.subset(coltbl, 1:2)) == (X = X[1:2], Y = Y[1:2], Z = Z[1:2])

    @test replacetable(coltbl, baz).tbl == baz
    @test CausalTables.replace(rowtbl; treatment = :X).treatment == :X

    # Errors
    @test_throws ArgumentError CausalTable(foo1, :X, :X, [:Z])
    @test_throws ArgumentError CausalTable(foo1, :X, :Z, [:Z])
    @test_throws ArgumentError CausalTable(foo1, :Z, :X, [:Z])

end

@testset "DataGeneratingProcess, no graphs" begin
    distseq = [
        :L1 => (; O...) -> Beta(1, 1),
        :L2 => (; O...) -> Multinomial(length(O[:L1]), O[:L1] / sum(O[:L1])),
        :A => (; O...) -> (@. Normal(O[:L1], 1)),
        :Y => (; O...) -> (@. Normal(O[:A] + 0.2 * O[:L2], 1))
    ]

    dgp = DataGeneratingProcess(distseq);
    foo = rand(dgp, 10)
    
    @test typeof(foo) == CausalTable
    @test Tables.columnnames(foo.tbl) == (:L1, :L2, :A, :Y)

    bar = condensity(dgp, foo, :A)
    baz = conmean(dgp, foo, :Y)
    @test nrow(foo.tbl) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}
    @test typeof(baz) <: Vector{T} where {T <: Real}
    @test baz == Tables.getcolumn(foo, :A) .+ 0.2 .* Tables.getcolumn(foo, :L2)
end

@testset "DataGeneratingProcess using dgp macro, no graphs" begin
    distseq = @dgp(
        L1 ~ Beta(1,1),
        L2 ~ Multinomial(length(:L1), :L1 / sum(:L1)),
        A ~ (@. Normal(:L1, 1)),
        Y ~ (@. Normal(:A + 0.2 * :L2, 1))
    )

    dgp = DataGeneratingProcess(distseq);
    foo = rand(dgp, 10)
    
    @test typeof(foo) == CausalTable
    @test Tables.columnnames(foo.tbl) == (:L1, :L2, :A, :Y)

    bar = condensity(dgp, foo, :A)
    baz = conmean(dgp, foo, :Y)
    @test nrow(foo.tbl) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}
    @test typeof(baz) <: Vector{T} where {T <: Real}
    @test baz == Tables.getcolumn(foo, :A) .+ 0.2 .* Tables.getcolumn(foo, :L2)

end


@testset "DataGeneratingProcess with graphs" begin

    distseq = Vector{Pair{Symbol, CausalTables.ValidDGPTypes}}([
        :L1 => (; O...) -> DiscreteUniform(1, 5),
        :L1_s => Sum(:L1),
        :A => (; O...) -> (@. Normal(O[:L1] + O[:L1_s], 1)),
        :A_s => Sum(:A),
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
    

    # Test the graph subsetting capabilities of CausalTable
    indices = 1:10
    baz = Tables.subset(foo, indices)
    @test baz.tbl == Tables.subset(foo.tbl, indices)
    @test nv(baz.graph) == nv(foo.graph[indices])

    # Test that duplicate indices throw an error
    @test_throws ArgumentError Tables.subset(foo, [2, 2])
end

@testset "DataGeneratingProcess with graphs using dgp macro" begin
    distseq = @dgp(
        L1 ~ DiscreteUniform(1, 5),
        L1_s = Sum(:L1, include_self = false),
        A ~ (@. Normal(:L1 + :L1_s, 1)),
        A_s = Sum(:A, include_self = false),
        Y ~ (@. Normal(:A + :A_s + 0.2 * :L1 + 0.05 * :L1_s, 1))
    )

    dgp = DataGeneratingProcess(n -> erdos_renyi(n, 3/n), distseq; controls = [:L1, :L1_s]);
    foo = rand(dgp, 100)
    @test typeof(foo) == CausalTable
    @test Tables.columnnames(foo.tbl) == (:L1, :L1_s, :A, :A_s, :Y)

    bar = condensity(dgp, foo, :A_s)
    foo_sum = summarize(foo)
    @test nrow(foo.tbl) == length(bar)
    @test typeof(bar) <: Vector{T} where {T <: UnivariateDistribution}   
    @test Tables.getcolumn(foo, :L1_s) == Tables.getcolumn(summarize(foo), :L1_s)
    

    # Test the graph subsetting capabilities of CausalTable
    indices = 1:10
    baz = Tables.subset(foo, indices)
    @test baz.tbl == Tables.subset(foo.tbl, indices)
    @test nv(baz.graph) == nv(foo.graph[indices])
end

@testset "Summary objects" begin
    @test get_var_to_summarize(Sum(:A)) == :A

    Random.seed!(1)
    distseq = @dgp(
        A ~ (@. Normal(0, 1)),
        A_sum = Sum(:A, include_self = false),
        A_sum2 = Sum(:A, include_self = true),
        A_max = Maximum(:A, include_self = false),
        A_max2 = Maximum(:A, include_self = true, use_inneighbors = false),
        A_min = Minimum(:A, include_self = false, use_inneighbors = false),
        A_min2 = Minimum(:A, include_self = true),
        A_prod = Product(:A, include_self = false),
        A_prod2 = Product(:A, include_self = true),
        F1 = Friends(),
        F2 = Friends(use_inneighbors = false),
        B ~ DiscreteUniform(1, 4),
        B_mode = Mode(:B, include_self = false),
        B_mode2 = Mode(:B, include_self = true),
        B_prod = Product(:B, include_self = false)
    )

    dgp = DataGeneratingProcess(n -> random_regular_graph(n, 5), distseq);
    data = rand(dgp, 10)
    data2 = CausalTables.summarize(data)
    tbl2 = CausalTables.summarize(data; keep_original = false)

    @test gettable(data) == gettable(data2)
    @test TableOperations.select(data2, :A_sum, :A_sum2, :A_max, :A_max2, :A_min, :A_min2, :A_prod, :A_prod2, :F1, :F2, :B_mode, :B_mode2, :B_prod) |> Tables.columntable == tbl2

    i = 1
    f = neighbors(data.graph, i)
    A = Tables.getcolumn(data, :A)
    A_samp = A[f]
    B = Tables.getcolumn(data, :B)
    B_samp = B[f]

    @test Tables.getcolumn(data, :F1)[i] == 5
    @test Tables.getcolumn(data, :F2)[i] == 5
    @test Tables.getcolumn(data, :A_sum)[i] == sum(A_samp)
    @test Tables.getcolumn(data, :A_sum2)[i] == sum(A_samp) + A[i]
    @test Tables.getcolumn(data, :A_prod)[i] ≈ prod(A_samp)
    @test Tables.getcolumn(data, :A_prod2)[i] ≈ prod(A_samp) * A[i]
    @test Tables.getcolumn(data, :B_mode)[i] == mode(B_samp)
    @test Tables.getcolumn(data, :B_mode2)[i] == mode([B_samp; B[i]])
    @test Tables.getcolumn(data, :A_max)[i] == maximum(A_samp)
    @test Tables.getcolumn(data, :A_max2)[i] == maximum([A_samp; A[i]])
    @test Tables.getcolumn(data, :A_min)[i] == minimum(A_samp)
    @test Tables.getcolumn(data, :A_min2)[i] == minimum([A_samp; A[i]])
    Tables.getcolumn(data, :B_prod)[i]
    @test Tables.getcolumn(data, :B_prod)[i] ≈ prod(B_samp)
end


@testset "DGP Exception throwing" begin

    ### Test that the DGP macro throws an error when it should ###
    
    # Test equality operator
    @test_throws ArgumentError CausalTables._parse_tilde(:(A = Normal(0, 1)))

    # Test LHS
    @test_throws ArgumentError CausalTables._parse_tilde(:(a() = Normal(0, 1))) 
    @test_throws ArgumentError CausalTables._parse_tilde(:(1 ~ Normal(0, 1)))
    @test_throws ArgumentError CausalTables._parse_tilde(:(1a ~ Normal(0, 1)))
    @test_throws ArgumentError CausalTables._parse_tilde(:(A; ~ Normal(0, 1)))
    @test_throws ArgumentError CausalTables._parse_tilde(:([]p23[p4] ~ Normal(0, 1)))
    
    ### Test the DGP constructor

    distseq = @dgp(
        L1 ~ DiscreteUniform(1, 5),
        L1_s = Sum(:L1, include_self = false),
        L1_s2 = Product(:L1, include_self = false),
        A ~ (@. Normal(:L1 + :L1_s, 1)),
        A_s = Sum(:A, include_self = false),
        Y ~ (@. Normal(:A + :A_s + 0.2 * :L1 + 0.05 * :L1_s, 1))
    )

    @test_throws ArgumentError DataGeneratingProcess(distseq; treatment = :L1, controls = [:L1])
    @test_throws ArgumentError DataGeneratingProcess(distseq; response = :L1, controls = [:L1])
    @test_throws ArgumentError DataGeneratingProcess(distseq; treatment = :X, controls = [:Y])

    # Test the RHS
    # TODO: Currently errors in DGP construct are deferred until rand or condensity is called.
    # Can we catch them earlier?
    bad = DataGeneratingProcess(@dgp(L ~ asdjfk))
    @test_throws ErrorException rand(bad, 10)
    tbl = CausalTable((L = [1, 2, 3],))
    @test_throws ErrorException condensity(bad, tbl, :L)

    bad1 = DataGeneratingProcess(distseq; treatment = :A, response = :Y, controls = [:L1])
    @test_throws ArgumentError rand(bad1, 10)
    bad2 = DataGeneratingProcess(n -> Graphs.path_graph(n), distseq, controls = [:L1_s2])
    badtab = rand(bad2, 5)
    @test_throws ErrorException condensity(bad2, badtab, :L1_s)
    @test_throws ErrorException condensity(bad2, badtab, :L1_s2)

    bad3 = DataGeneratingProcess(n -> throw(ErrorException("bad function")), distseq)
    @test_throws ErrorException rand(bad3, 10)
end
