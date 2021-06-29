using Test, JuliaDispatch

include("../src/_aux.jl")
include("../src/_aux2.jl")

@testset "Test read aux" begin
    auxDict_python = aux(file="data/00001.aux")
    auxDict_julia = aux2(file="data/00001.aux")

    for (k, v) in auxDict_julia
        @test k in keys(auxDict_python)
        @test v == auxDict_python[k]
    end
end
