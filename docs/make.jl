using UnitfulTensors, Documenter, DocumenterMermaid

DocMeta.setdocmeta!(UnitfulTensors, :DocTestSetup, :(using UnitfulTensors); recursive=true)

makedocs(sitename="UnitfulTensors.jl",
         modules=[UnitfulTensors, UnitfulTensors.FastQuantities],
         remotes=nothing,
         pages = ["Home" => "index.md", "guide.md", "api.md"],
         checkdocs=:exports)
         
deploydocs(
    repo = "github.com/anonymous-shrew/UnitfulTensors.jl.git",
)
