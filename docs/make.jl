using Documenter, CausalTables

makedocs(
    sitename="CausalTables.jl",
    pages = [
        "Home" => "index.md",
        "Turning data into a `CausalTable`" =>              "man/formatting.md",
        "Generating data for statistical experiments" =>    "man/generating-data.md",
        "Approximating ground truth causal estimands" =>    "man/estimands.md",
        "Computing ground truth conditional distributions" =>    "man/ground-truth.md",
        "Network summaries" =>                 "man/network-summaries.md",
        "API" =>                                "man/api.md",
    ]
)

deploydocs(
    repo = "github.com/salbalkus/CausalTables.jl.git"
)