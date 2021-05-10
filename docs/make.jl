using ShellModel
using Documenter

DocMeta.setdocmeta!(ShellModel, :DocTestSetup, :(using ShellModel); recursive=true)

makedocs(;
    modules=[ShellModel],
    authors="SotaYoshida <s.yoshida@nt.phys.s.u-tokyo.ac.jp> and contributors",
    repo="https://github.com/SotaYoshita/ShellModel.jl/blob/{commit}{path}#{line}",
    sitename="ShellModel.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://SotaYoshita.github.io/ShellModel.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SotaYoshita/ShellModel.jl.git",
)
