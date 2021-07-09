using JuliaDispatch
using Test

@testset "JuliaDispatch.jl" begin
    include("TestExpressionParser.jl")
    include("TestSnapshot.jl")
end
