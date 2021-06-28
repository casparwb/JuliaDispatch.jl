using Conda, PyCall, Test

include("_dispatch.jl")
cd("experiments\\shu-osher\\julia")
function test_patch()
    """ Test patch from shu-osher snapshot 1 """

    top = "C:\\Users\\caspa\\Documents\\DIspatch\\dispatch\\"

    pushfirst!(PyVector(pyimport("sys")."path"), top*"utilities\\python");

    dispatch = pyimport("dispatch")

    data = "../data"
    run = ""

    i = 0

    snap_julia = snapshot(i)["patches"][1]
    snap_python = dispatch.snapshot(i,run=run,data=data, verbose=0).patches[0]
    attrs = collect(keys(snap_julia))
    for kpy in keys(snap_python)
        if !(String(kpy) in attrs) && !startswith(String(kpy), "__")
            println(kpy)
        end
    end


end

function test_snapshot()

    """ Test patch from shu-osher snapshot 1 """

    top = "C:\\Users\\caspa\\Documents\\DISPATCH_julia"

    pushfirst!(PyVector(pyimport("sys")."path"), top*"utilities\\python");

    dispatch = pyimport("dispatch")

    data = "../data"
    run = ""

    i = 0

    snap_julia = snapshot(i)
    snap_python = dispatch.snapshot(i,run=run,data=data, verbose=0)
    attrs = sort(collect(keys(snap_julia)))
    for kpy in sort(collect(keys(snap_python)))
        if !(String(kpy) in attrs) && !startswith(String(kpy), "__")
            println(kpy)
        end
    end
