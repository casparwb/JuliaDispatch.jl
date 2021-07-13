using Pkg
Pkg.activate(".")
using Documenter
using JuliaDispatch

makedocs(
    modules=[JuliaDispatch, JuliaDispatch.Dispatch, JuliaDispatch.Graphics, 
             JuliaDispatch.Select, JuliaDispatch.Buffers, JuliaDispatch.Utils],
    sitename = "JuliaDispatch.jl",
    authors = "Caspar William Bruenech",
    format = Documenter.HTML(),
    build= "build",
    pages = ["JuliaDispatch: Analysis and Visualization Tools for Dispatch" => "index.md",
            "Getting Started" => "installation.md",
            "Quick Start" => "quickstart.md",
            "Dispatch" => "dispatch.md",
            "Graphics" => "graphics.md",
            "Select" => "select.md",
            "Buffers" => "buffers.md",
            "Utils" => "utils.md"]
        )

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
# deploydocs(
#     repo = "github.com/casparwb/JuliaDispatch.jl.git"
# )

deploydocs(;
    repo="github.com/casparwb/JuliaDispatch.jl.git"
)
