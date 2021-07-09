module TestExpressionParser

    using JuliaDispatch.Dispatch, Test

    println("Testing expression parser (_expr.jl)")

    isfile("data/orz/data/params.jld") && rm("data/orz/data/params.jld")
    isfile("data/orz/data/00100/snapshot.jld") && rm("data/orz/data/00100/snapshot.jld")
    isfile("data/orz/data/00100/patches.jld") && rm("data/orz/data/00100/patches.jld")

    snap = snapshot(100, data="data/orz/data")

    expression = "π/4*sqrt(bx^2 + by^2 + bz^2)"

    @testset "Expression parser test" begin
        @testset "patch $patch_id" for patch_id = 1:length(snap["patches"])
            patch = snap["patches"][patch_id]
            var_parsed = patch["var"](expression)
            var_manual = π/4*sqrt.( patch["var"]("bx") .^ 2 .+ patch["var"]("by") .^ 2 .+ patch["var"]("bz") .^ 2 )
            @test var_parsed == var_manual
        end
    end
end