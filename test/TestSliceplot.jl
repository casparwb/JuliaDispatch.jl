module TestGraphics
    using JuliaDispatch.Dispatch, JuliaDispatch.Graphics

    snap = snapshot(0, data="data/orz/data");
    sliceplot(snap, z=0.1, iv="d");

end