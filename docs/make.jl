using FVUnraveling
using Documenter

DocMeta.setdocmeta!(FVUnraveling, :DocTestSetup, :(using FVUnraveling); recursive=true)

makedocs(;
    modules=[FVUnraveling],
    authors="Malte Krug <malte.krug@uni-ulm.de> and contributors",
    sitename="FVUnraveling.jl",
    format=Documenter.HTML(;
        canonical="https://krugx.github.io/FVUnraveling.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/krugx/FVUnraveling.jl",
    devbranch="main",
)
