using Test
using CausalTables
using Tables
using TableOperations
import StatsBase: mode
using DataFrames
using Graphs
using Distributions
using Random


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
    baz = (X = [4, 5], Y = ["foo", "bar"], Z = [0.1, 0.2])
    gettable(CausalTables.replace(rowtbl; tbl = baz))
    @test gettable(CausalTables.replace(rowtbl; tbl = baz)) == baz
    @test gettable(Tables.subset(coltbl, 1:2)) == (X = X[1:2], Y = Y[1:2], Z = Z[1:2])

    # Errors
    @test_throws ArgumentError CausalTable(foo1, :X, :X, [:Z])
    @test_throws ArgumentError CausalTable(foo1, :X, :Z, [:Z])
    @test_throws ArgumentError CausalTable(foo1, :Z, :X, [:Z])

    tbl4 = CausalTable(foo1, :W, :Y, [:Z])
    @test_throws ErrorException("Treatment variable not contained in the data. Note: If response is a summary over a network (contained within tbl.summaries), make sure that you call `summary(tbl::CausalTable)` on your table before calling `gettreatment`.") gettreatment(tbl4)
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

@testset "DataGeneratingProcess using dgp macro, no graphs" begin
    distseq = @dgp(
        L1 ~ DiscreteUniform(1, 5),
        A ~ (@. Normal(:L1, 1)),
        Y ~ (@. Normal(:A + 0.2 * :L1, 1))
    )

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
end

@testset "DataGeneratingProcess with graphs using dgp macro" begin
    distseq = @dgp(
        L1 ~ DiscreteUniform(1, 5),
        L1_s = Sum(:L1),
        A ~ (@. Normal(:L1 + :L1_s, 1)),
        A_s = Sum(:A),
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

@testset "test all summary functions" begin
    Random.seed!(1)
    distseq = @dgp(
        A ~ (@. Normal(0, 1)),
        A_sum = Sum(:A, include_self = false),
        A_sum2 = Sum(:A, include_self = true),
        A_max = Maximum(:A, include_self = false),
        A_max2 = Maximum(:A, include_self = true),
        A_min = Minimum(:A, include_self = false),
        A_min2 = Minimum(:A, include_self = true),
        A_prod = Product(:A, include_self = false),
        A_prod2 = Product(:A, include_self = true),
        F = Friends(),
        B ~ Binomial(4, 0.5),
        B_mode = Mode(:B, include_self = false),
        B_mode2 = Mode(:B, include_self = true)
    )

    dgp = DataGeneratingProcess(n -> random_regular_graph(n, 5), distseq);
    data = rand(dgp, 10)
    data2 = CausalTables.summarize(data)
    tbl2 = CausalTables.summarize(data; keep_original = false)

    @test gettable(data) == gettable(data2)
    @test TableOperations.select(data2, :A_sum, :A_sum2, :A_max, :A_max2, :A_min, :A_min2, :A_prod, :A_prod2, :F, :B_mode, :B_mode2) |> Tables.columntable == tbl2

    i = 1
    f = neighbors(data.graph, i)
    A = Tables.getcolumn(data, :A)
    A[i]
    A_samp = A[f]
    B_samp = Tables.getcolumn(data, :B)[f]
    @test Tables.getcolumn(data, :F)[i] == 5
    @test Tables.getcolumn(data, :A_sum)[i] == sum(A_samp)
    @test Tables.getcolumn(data, :A_sum2)[i] == sum(A_samp) + A[i]
    @test Tables.getcolumn(data, :A_prod)[i] == prod(A_samp)
    @test Tables.getcolumn(data, :A_prod2)[i] ≈ prod(A_samp) * A[i]
    @test Tables.getcolumn(data, :B_mode)[i] == mode(B_samp)
    @test Tables.getcolumn(data, :B_mode2)[i] == mode(B_samp)
    @test Tables.getcolumn(data, :A_max)[i] == maximum(A_samp)
    @test Tables.getcolumn(data, :A_max2)[i] == maximum([A_samp; A[i]])
    @test Tables.getcolumn(data, :A_min)[i] == minimum(A_samp)
    @test Tables.getcolumn(data, :A_min2)[i] == minimum([A_samp; A[i]])
end


