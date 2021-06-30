module TestModules
    using Test
    @testset "test import" begin
        @test using JuliaDispatch.Dispatch
        @test using JuliaDispatch.Select
        @test using JuliaDispatch.Graphics
        @test using JuliaDispatch.Buffers
        @test using JuliaDispatch.Analysis
        @test using JuliaDispatch.Utils
    end

end