module TestSnapshot

    using JuliaDispatch.Dispatch
    using PyCall, Test
    

    snap0 = snapshot(0, data="data/orz/data", suppress=true, progress=false)     # load first snapshot 
    snap100 = snapshot(100, data="data/orz/data", suppress=true, progress=false) # load 100th snapshot
    
    # import the data from the Python analysis
    np = pyimport("numpy")
    snap0py = PyDict(np.load("data/snap0.npy", allow_pickle=true)[1])
    snap100py = PyDict(np.load("data/snap100.npy", allow_pickle=true)[1])

    # test that snapshots from both languages have the same properties (keys) 
    @testset "Snapshot keys $key0" for (key0, key100) in zip(collect(snap0py["keys"]), collect(snap100py["keys"]))
            @test key0 in keys(snap0)
            @test key100 in keys(snap100)
    end

    # test that the number of patches is the same
    @testset "Number of patches" begin
        npatches0 = length(snap0["patches"])
        npatches100 = length(snap100["patches"])
        
        @test npatches0 == snap0py["npatches"]
        @test npatches100 == snap100py["npatches"]
    end

    # compare a number of quantities for a selection of patches
    @testset verbose=true "ivs" begin
        patches = [patch for patch in keys(snap0py) if startswith(patch, "patch")]
        @testset "patch $patch" for patch in patches
            ivs = keys(snap0py[patch])
            patchID = parse(Int, split(patch, "_")[2])
            patch0 = findall(x -> x["id"] == patchID, snap0["patches"])[1]
            patch100 = findall(x -> x["id"] == patchID, snap100["patches"])[1]
            @testset "iv $iv" for iv in ivs
                iv_python_0 = snap0py[patch][iv] |> copy
                iv_python_100 = snap100py[patch][iv] |> copy

                iv_julia_0 = snap0["patches"][patch0]["var"](iv) |> copy
                iv_julia_100 = snap100["patches"][patch100]["var"](iv) |> copy

                replace!(iv_python_0, NaN => 0.0)
                replace!(iv_python_100, NaN => 0.0)
                replace!(iv_julia_0, NaN => 0.0)
                replace!(iv_julia_0, NaN => 0.0)
                
                @testset "Snapshot 0" begin
                    @test iv_python_0 ≈ iv_julia_0
                end # testset

                @testset "Snapshot 100" begin
                    @test iv_python_100 ≈ iv_julia_100
                end # testset
            end # for 
        end # for 
    end # testset 

    # test auxiliary parameters
    auxvars = snap100["aux"]["select"]
    @testset "Aux quantity $auxIV" for auxIV in auxvars
        @test auxIV in snap100py["aux"]["select"]
    end # testset
    
    patches = [patch for patch in keys(snap100py) if startswith(patch, "patch")]
    @testset "patch $patch" for patch in patches
        patchID = parse(Int, split(patch, "_")[2])
        patch_idx = findall(x -> x["id"] == patchID, snap100["patches"])[1]
        patch100_j = snap100["patches"][patch_idx]
        patch100_p = snap100py[patch]

        @test patch100_j["aux"]["id"] == convert(Int, patch100_p["aux"]["id"])
        @test patch100_j["aux"]["filename"] == patch100_p["aux"]["filename"]
        @test patch100_j["aux"]["version"] == convert(Int, patch100_p["aux"]["version"])

        @testset "aux var $var" for var in keys(patch100_p["aux"]["vars"])
            auxvar_j = patch100_j["aux"]["vars"][var]
            auxvar_p = patch100_p["aux"]["vars"][var]

            auxvar_j_v = patch100_j["var"](var)
            auxvar_p_v = patch100_p["aux"]["vars"][var]["v"] |> copy

            # replace NaNs, as NaN !== NaN
            replace!(auxvar_j_v, NaN => 0.0)
            replace!(auxvar_p_v, NaN => 0.0)
            
            @test auxvar_j["rank"] == auxvar_p["rank"]
            @test auxvar_j["name"] == auxvar_p["name"]
            @test all(auxvar_j["shape"] .== auxvar_p["shape"])
            @test auxvar_j["type"] == auxvar_p["type"]
            @test auxvar_j_v ≈ auxvar_p_v
        end #testset
    end # testset


    # # test patch axes
    # patches = [patch for patch in keys(snap0py) if startswith(patch, "patch")]
    # @testset "patch $patch" for patch in patches
    #     ivs = keys(snap0py[patch])
    #     patchID = parse(Int, split(patch, "_")[2])
    #     patch0 = findall(x -> x["id"] == patchID, snap0["patches"])[1]
    #     patch100 = findall(x -> x["id"] == patchID, snap100["patches"])[1]


end # module