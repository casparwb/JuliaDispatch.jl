module TestExpressionParser

    using JuliaDispatch.Dispatch, Test

    @info "Testing expression parser (_expr.jl)"

    snap = snapshot(100, data="data/orz/data", suppress=true, progress=false)

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