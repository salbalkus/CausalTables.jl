using Documenter, CausalTables

makedocs(
    sitename="CausalTables.jl",
    pages = [
        "Home" => "index.md".
        "Getting Started" => "getting_started.md",
        "Turning data into a `CausalTable`" => "formatting.md",
        "Generating Data for Statistical Experiments" => "generating_data.md",
        "Computing Ground Truth of Causal Parameters" => "ground_truth.md"
        "User-Defined Network Summaries" => "network_summaries.md",
    ]
)

deploydocs(
    repo = "github.com/salbalkus/CausalTables.jl.git"
)