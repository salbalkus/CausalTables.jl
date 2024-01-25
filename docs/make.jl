using Documenter, CausalTables

makedocs(
    sitename="CausalTables.jl",
    pages = [
        "Home" => "index.md",
        "Getting Started" =>                                "man/getting-started.md",
        "Turning data into a `CausalTable`" =>              "man/formatting.md",
        "Generating Data for Statistical Experiments" =>    "man/generating-data.md",
        "Computing Ground Truth of Causal Parameters" =>    "man/ground-truth.md",
        "Network Summaries" =>                 "man/network-summaries.md",
        "API" =>                                "man/api.md",
    ]
)

deploydocs(
    repo = "github.com/salbalkus/CausalTables.jl.git"
)