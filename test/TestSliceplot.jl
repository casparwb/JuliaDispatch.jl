module TestGraphics
    using Test
    using JuliaDispatch.Dispatch, JuliaDispatch.Graphics

    snap = snapshot(0, data="data/orz/")
    @test sliceplot(snap, z=0.1, iv="d")

end