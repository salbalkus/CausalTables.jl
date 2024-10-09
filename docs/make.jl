using Documenter, CausalTables

makedocs(
    sitename="CausalTables.jl",
    pages = [
        "Home" => "index.md",
        "Turning data into a `CausalTable`" =>              "man/formatting.md",
        "Generating Data for Statistical Experiments" =>    "man/generating-data.md",
        "Computing Ground Truth Conditional Distributions" =>    "man/ground-truth.md",
        "Approximating Ground Truth Causal Estimands" =>    "man/estimands.md",
        "Network Summaries" =>                 "man/network-summaries.md",
        "API" =>                                "man/api.md",
    ]
)

deploydocs(
    repo = "github.com/salbalkus/CausalTables.jl.git"
)