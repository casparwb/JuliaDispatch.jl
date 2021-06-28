using Revise
include("_dispatch.jl")
include("graphics/_graphics.jl")
include("select/_select.jl")

get_patch(i=1) = snapshot()["patches"][i];
