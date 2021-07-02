module TestSnapshot

    using JuliaDispatch.Dispatch
    using PyCall, Test, JLD
    
    isfile("data/orz/data/params.jld") && rm("data/orz/data/params.jld")
    isfile("data/orz/data/00000/snapshot.jld") && rm("data/orz/data/00000/snapshot.jld")
    isfile("data/orz/data/00000/patches.jld") && rm("data/orz/data/00000/patches.jld")
    isfile("data/orz/data/00100/snapshot.jld") && rm("data/orz/data/00100/snapshot.jld")
    isfile("data/orz/data/00100/patches.jld") && rm("data/orz/data/00100/patches.jld")

    snap0 = snapshot(0, data="data/orz/data")
    rm("data/orz/data/params.jld")
    snap100 = snapshot(100, data="data/orz/data")
    
    np = pyimport("numpy")
    snap0data = PyDict(np.load("data/snap0.npy", allow_pickle=true)[1])
    snap100data = PyDict(np.load("data/snap100.npy", allow_pickle=true)[1])

    @testset "Snapshot keys" begin
        for (key0, key100) in zip(keys(snap0), keys(snap100))
            @test key0 in snap0data["keys"]
            @test key100 in snap100data["keys"]
        end
    end

    @testset "Number of patches" begin
        npatches0 = length(snap0["patches"])
        npatches100 = length(snap100["patches"])
        
        @test npatches0 == snap0data["npatches"]
        @test npatches100 == snap100data["npatches"]
    end

    @testset verbose=true "ivs" begin
        patches = [patch for patch in keys(snap0data) if startswith(patch, "patch")]
        @testset "patches $patch" for patch in patches
            ivs = keys(snap0data[patch])
            patchID = parse(Int, split(patch, "_")[2])
            patch0 = findall(x -> x["id"] == patchID, snap0["patches"])[1]
            patch100 = findall(x -> x["id"] == patchID, snap100["patches"])[1]
            @testset "patch ivs $iv" for iv in ivs
                iv_python_0 = snap0data[patch][iv]
                iv_python_100 = snap100data[patch][iv]

                iv_julia_0 = snap0["patches"][patch0]["var"](iv)
                iv_julia_100 = snap100["patches"][patch100]["var"](iv)

                @test iv_python_0 ≈ iv_julia_0
                @test iv_python_100 ≈ iv_julia_100
            end
        end
    end
end