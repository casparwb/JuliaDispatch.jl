using PyCall, Printf, Mmap, StaticArrays, JLD, Unitful, ProgressBars
using FStrings

using JuliaDispatch.Utils

include("_aux2.jl")
include("_dispatch_grid.jl")
include("dispatch_utils.jl")
include("../expr/_expr.jl")

"""
    snapshot(iotu::Int; run::String, data::String, verbose::Int)

Parses patches and produces a dictionary with all properties of snapshot
`iout`.

#Arguments:
- `iout::Int`, snapshot ID

#Kwargs:
- `data::String`, path to data (snapshots), default `"..\\data"`
- `run::String`, path to snapshots folders relative to data, default `""`

#Returns:
- Dictionary of snapshot
"""
function snapshot(iout=0; run="", data="../data", verbose = 0, copy = false, memmap = 1)

    ### wrapper for f90nml reader ###
    #read_nml(file) = pyimport("f90nml").read(file)

    ### find data- and rundirs ###
    # rundir = _dir(data, run)
    statedict = Dict{String, Any}()
    rundir = _dir(data, run)
    if !isdir(rundir)
        println("directory $(rundir) does not exist")
        return nothing
    end

    statedict["copy"] = copy

    #datadir = _dir(rundir, "$(@sprintf("%05d", iout))")
    datadir = joinpath(rundir, "$(@sprintf("%05d", iout))")

    if !isdir(datadir)
        println("directory $(datadir) does not exist")
        return nothing
    end

    #=
    statedict (temporary name) is the equivalent of the snapshot class in
    the python implementation
    =#
    statedict["datadir"] = datadir
    statedict["rundir"] = rundir

    ### add properties from namelist ###
    file = _file(rundir, "params.nml")

    #params_list = read_nml(file)
    statedict["params_list"] = read_nml(file)
    _add_nml_list_to(statedict, statedict["params_list"])

    ### finding snapshot files ###
    files = [f for f in readdir(datadir) if endswith(f, "snapshot.nml")]

    println("Number of snapshot.nml files in $datadir: $(length(files))")
    if length(files) == 0
        println("Directory $datadir has no snapshot.nml file")
        return nothing
    end

    file_dict = []
    @inbounds for f in files
        file = _file(datadir, f)
        nml_list = read_nml(file)
        statedict["nml_list"] = nml_list
        _add_nml_list_to(statedict, nml_list)

        if "idx_nml" in keys(nml_list)
            idx_dict = nml_list["idx_nml"]
            for (k, v) in idx_dict
                idx_dict[k] += 1
            end

            idx = Dict()
            idx["dict"] = idx_dict
            idx["vars"] = Dict()
            for (k, v) in idx["dict"]
                if !(v in keys(idx["vars"]))
                    v > 0 ? idx["vars"][v] = k : nothing
                end
            end

            statedict["idx"] = idx
            for (k, v) in idx["dict"]
                statedict["idx"][k] = v
            end

            statedict["keys"] = []
            for (k, v) in statedict["idx"]["vars"]
                push!(statedict["keys"], k)
                push!(statedict["keys"], v)
            end
        end
    end



    ### Check if snapshots.dat ###
    verbose > 1 && println("checking if there is a snapshots.dat")
    try
        statedict["datafile"] = _file(rundir, "snapshots.dat")
        statedict["datafiled"] = open(statedict["datafile"], "rb")
    catch
        nothing
    end
    statedict["dict"] = statedict["nml_list"]["snapshot_nml"]

    ### dress snapshot with attributes for snapshot_nml ###

    """ add patches as a list of dicts """
    verbose > 1 && println("add patches as a list of dicts") 

    statedict["patches"] = Dict[]#Array{Dict, 1}([])
    files = [f for f in readdir(datadir) if endswith(f, "_patches.nml")]


    if !haskey(statedict, "mpi_size")
        statedict["mpi_size"] = 1
    end

    patch_dict = read_patch_metadata(iout, run, data, statedict["mpi_size"], verbose=verbose)

    if isnothing(patch_dict)
        println("could not find patch metadata")
        return nothing
    end

    ids = sort(collect(keys(patch_dict)))
    for id in ProgressBar(1:length(ids))
        p = _patch2(id, patch_dict["$id"], statedict)
        # p = _patch2(parse(Int, id), patch_dict[id], statedict)
        _add_axes(statedict, p)
        push!(statedict["patches"], p)
        
        if verbose == 2 && haskey(p, "idx")
            data = p["var"]("d")
            dmax = maximum(data)
            println("id: $(p["id"])   pos: $(p["position"]) $dmax")
        elseif verbose == 3
            println("id: $(p["id"])")
            for iv in range(p["nv"])
                data = p["var"](iv)
                vmin = parse(Float64, minimum(data))
                vmax = parse(Float64, maximum(data))
                var = p["idx"]["vars"][iv]
                println(@sprintf("%d :  min = %10.3e    max = %10.3e", var, vmin, vmax))
            end # do
        elseif verbose == 4
            attributes(p)
        end # if
    end # for

    verbose == 1 && println("    added $(length(statedict[string(:patches)])) patches")

    return statedict

end

""" 
    attributes(patch)

Pretty-print patch attributes
"""
function attributes(patch)
    id = patch["id"]
    println("id: $id")

    for (k, v) in patch
        if k != "data"
            #print(f"{k:>12}: {v}")
        end # if 
    end # for
end # function

""" Add axis values to patches """
function _add_axes(snap, patch)

    first = patch["llc_cart"] .- patch["ng"].*patch["ds"]

    n = patch["gn"]
    ng = patch["ng"]

    if patch["no_mans_land"]
        first .= first .+ 0.5*patch["ds"]
    end
    patch["x"] = first[1] .+ patch["ds"][1]*collect(0:n[1]-1)
    patch["y"] = first[2] .+ patch["ds"][2]*collect(0:n[2]-1)
    patch["z"] = first[3] .+ patch["ds"][3]*collect(0:n[3]-1)


    patch["xi"] = patch["x"][(ng[1] + 1):end-ng[1]]
    patch["yi"] = patch["y"][(ng[2] + 1):end-ng[2]]
    patch["zi"] = patch["z"][(ng[3] + 1):end-ng[3]]
    patch["xs"] = patch["x"] .- 0.5*patch["ds"][1]
    patch["ys"] = patch["y"] .- 0.5*patch["ds"][2]
    patch["zs"] = patch["z"] .- 0.5*patch["ds"][3]
    patch["xyz"] = [patch["x"], patch["y"], patch["z"]]
    patch["xyzi"] = [patch["xi"], patch["yi"], patch["zi"]]

    # add geometric factors to patch (important for spherical/cylindrical coords.)
    if snap["mesh_type"] == 1
        patch["geometric_factors"] = nothing
    else
        patch["geometric_factors"] = GeometricFactors(patch)
    end
end

"""
    read_nml(file; verbose=0)

Read a cached namelist file if it exists, else create it.
"""
function read_nml(file; verbose=0)
    jldfile = replace(file, ".nml" => ".jld")

    if ispath(jldfile) && ctime(jldfile) > ctime(file)
        verbose > 0 && println("reading file $jldfile")
        nml_list = JLD.load(jldfile)
    else
        verbose > 0 && println("reading file $file")
        f90nml(file) = pyimport("f90nml").read(file)
        nml_list = f90nml(file)
        verbose > -1 && println("caching metadata to $jldfile")
        try
            JLD.save(jldfile, nml_list)
        catch
            nothing
        end # try
    end # if

    return nml_list
end # function


function _parse_namelist(items)

    pos = items[2:end]
    if length(pos) == 1
        pos[1] = replace(pos[1], "3*" => "")
        pos = [pos[1], pos[1], pos[1]]
    elseif length(pos) == 2
        if findfirst("2*", pos[1]) !== nothing
            pos[1] = replace(pos[1], "2*" => "")
            pos = [pos[1], pos[1], pos[2]]
        end
        if findfirst("2*", pos[2]) !== nothing
            pos[2] = replace(pos[2], "2*" => "")
            pos = [pos[1], pos[2], pos[2]]
        end
    end
    # choose right data type for the resulting array
    try
        pos = [parse(Int, p) for p in pos]
    catch
        pos = [parse(Float64, p) for p in pos]
    end

    return pos
end

function _dir(dir, subdir)

    p = endswith(dir, '/') ? dir*subdir : dir*'/'*subdir
    if endswith(p, '/')
        return p
    else
        return p*'/'
    end
end


function _file(dir, file)
    endswith(dir, "\\") ? p = dir*file : dir*"\\"*file

    p = joinpath(dir, file)

end

""" Add all namelists as dictionairy entries """
function _add_nml_list_to(dict::Dict, nml::Dict)

    for (key, nml_dict) in nml
        if key == "snapshot_nml"
            _add_nml_to(dict, nml_dict)
        else
            name = replace(
                   replace(
                   String(key), "_nml" =>""), "_params" => "")

            dict[name] = Dict()
            if typeof(nml_dict) <: AbstractArray
                println("WARNING: more than one "*key)
                nml_dict = nml_dict[1]
            end
            _add_nml_to(dict[name], nml_dict)

        end
    end
end

""" Add dictionairy entries in d to dict """
function _add_nml_to(dict::Dict, subdict::Dict)
    for (key, value) in subdict
        dict[key] = value
    end
end


function _patch2(id, patch_dict, snap; memmap=1, verbose=0)
    patch = Dict()
    patch["id"] = id
    patch["memmap"] = memmap

    # add general attributes from snapshot.nml
    for (k, v) in snap["dict"]
        patch[k] = v
    end 

    # add per-patch attributes from parsing
    for (k, v) in patch_dict
        patch[k] = v
    end 

    # if a short variant of the patch metadata does not has ds, compute it
    if !haskey(patch_dict, "ds")
        patch["ds"] = patch["size"] ./ patch["n"]
    end 

    if !patch["guard_zones"]
        patch["li"][:] .= 1
        patch["ui"][:] .= patch["n"]
    end 

    # add idx attribute
    if haskey(snap, "idx")
        patch["idx"] = snap["idx"]
        patch["idx"]["h"] = _h(patch)
    end

    if haskey(patch, "size") && haskey(patch, "n")
        if !haskey(patch, "ds")
            patch["dsx"] = patch["size"]/patch["n"]
        end # if
    end # if 

    if haskey(patch, "size") && haskey(patch, "position")#"size" in keys(patch) && "position" in keys(patch)
        llc = patch["position"] - patch["size"]/2.0
        urc = patch["position"] + patch["size"]/2.0

        patch["extent"] = reshape([llc[2] urc[2] llc[3] urc[3]
                                    llc[3] urc[3] llc[1] urc[1]
                                    llc[1] urc[1] llc[2] urc[2]], 3, 4)
        # patch["extent"] = SMatrix{3, 4, Float16}(llc[2], urc[2], llc[3], urc[3],
        #                                         llc[3], urc[3], llc[1], urc[1],
        #                                         llc[1], urc[1], llc[2], urc[2])
        patch["llc_cart"] = llc

    end

    if !haskey(patch, "units")
        if haskey(snap, "units")
            patch["units"] = snap["units"]
            if !haskey(patch, "u")
                patch["units"]["u"] = patch["units"]["l"]/patch["units"]["t"]
            end
            if !haskey(patch, "d")
                patch["units"]["d"] = patch["units"]["m"]/patch["units"]["l"]^3
            end
            if !haskey(patch, "p")
                patch["units"]["p"] = patch["units"]["d"]/patch["units"]["u"]^2
            end
            if !haskey(patch, "e")
                patch["units"]["e"] = patch["units"]["m"]/patch["units"]["u"]^2
            end
            if !haskey(patch, "b")
                patch["units"]["b"] = patch["units"]["u"]*sqrt(4π*patch["units"]["d"])
            end
        end
    end

    # modify `mesh_type` from integer to string for readability
    if haskey(patch, "mesh_type")
        if patch["mesh_type"] == 1
            patch["mesh_type"] = "Cartesian"
        elseif patch["mesh_type"] == 2
            patch["mesh_type"] = "spherical"
        elseif patch["mesh_type"] == 3
            patch["mesh_type"] = "cylindrical"
        end
    end
    if strip(snap["io"]["method"]) == "legacy"
        patch["filename"] = snap["rundir"]*
                        "$(@sprintf("%05d/%05d", patch["iout"], patch["id"])).dat"
        patch["var"] = _var(patch, patch["filename"], snap)
    elseif strip(snap["io"]["method"]) == "background"
        patch["filename"] = snap["rundir"]*"$(@sprintf("%05d/snapshot_%05d.dat", patch["iout"], patch["rank"]))"
        patch["var"] = _var(patch, patch["filename"], snap, verbose=verbose)
    else
        patch["var"] = _var(patch, snap["datafiled"], snap)
    end

    """ add a comprehensive set of variable keys """
    patch["aux"] = Dict()
    patch["data"] = Dict()
    patch["keys"] = Dict()
    patch["keys"]["letters"] = collect(keys(snap["idx"]["dict"]))
    patch["keys"]["numbers"] = collect(values(snap["idx"]["dict"]))
    patch["keys"]["known"] = ["d","lnd","logd","ux","uy","uz","u1","u2","u3",
                                "ee","E","T","eth", "Eth", "ekin", "emag", "S"]


    # attach an aux filename, if one exists for this task
    auxfile = "$(@sprintf("%05d", id)).aux"

    auxfile = _pathjoin(snap["datadir"], auxfile)
    patch["auxfile"] = auxfile

    # read the information from the aux file and add as attribute
    patch["aux"] = aux(id = id, rundir=snap["rundir"], datadir=snap["datadir"],
                        io = snap["iout"], file = auxfile, verbose=verbose)
    patch["keys"]["aux"] = collect(keys(patch["aux"]["vars"]))

    # for patch id = 1, add items with les than 100 elements to snapshot.aux
    if id == 1
        for (k, v) in patch["aux"]["vars"]
            if prod(v["shape"]) < 100
                snap["aux"]["name"] = v["v"]
            end # if 
        end # for
    end # if 

    """ Collect all keys in a single array """
    all = []
    for key_list in values(patch["keys"])
        append!(all, key_list)
    end
    patch["all_keys"] = all

    return patch
end # function



function _var(patch, filed, snap; verbose = 0, copy = nothing)

    bytes = Int64(4*prod(patch["ncell"]))

    shape = tuple(patch["ncell"]...)
    patch["offset"] = []

    # p["ip"] is the offset in the file; ranging from 0 to ntotal-1
    if "cartesian" in keys(snap)
        nrank = prod(snap["cartesian"]["mpi_dims"])
        ntask = prod(snap["cartesian"]["per_rank"])
    else
        nrank = 1
        ntask = patch["ntotal"]
    end
    patch["ip"] = (patch["id"] - 1 - patch["rank"]) ÷ nrank + patch["rank"]*ntask

    # verbose == 5 ? println("id: $(patch["id"])") : nothing

    for iv = 0:patch["nv"] - 1
        if strip(snap["io"]["method"]) == "background"
            offset = iv + (patch["record"]-1)*patch["nv"]
        else
            if patch["ioformat"] == 5
                offset = patch["ip"] + iv*patch["ntotal"]
                offset += patch["iout"]*patch["ntotal"]*patch["nv"]

            elseif (patch["ioformat"] >= 6 && patch["ioformat"] < 10) || patch["ioformat"] == 15
                offset = iv + patch["ip"]*patch["nv"]
                offset += patch["iout"]*patch["ntotal"]*patch["nv"]

            elseif patch["ioformat"] >= 10
                offset = patch["ip"] + iv*patch["ntotal"]
                offset += patch["iout"]*patch["ntotal"]*patch["nv"]
            else
                offset = iv
            end
        end
        offset *= bytes
        push!(patch["offset"], offset)
    end

    function mem(iv; verbose = 0)
        """
        Translate alphabetic variable keys to numeric
        """

        verbose == 1 && println("mem($iv)")

        if typeof(iv) == typeof("d")
            iv = patch["idx"][iv]
        end

        # iszero(iv) ? iv += 1 : nothing
        if patch["memmap"] == 1
            v = Mmap.mmap(filed, Array{Float32, length(shape)}, shape, patch["offset"][iv])
        end

    end # end mem

    function dnup(q, shift = 1, axis = 1)
        shift == 1 ? i = 1 : i = -1

        if axis == 1 && size(q)[1] > 1
            f = (q + circshift(q, (shift, 0, 0)))*0.5
            i == -1 ? f[end,:,:] = q[end,:,:] : f[i,:,:] = q[i,:,:]
        elseif axis == 2 && size(q)[2] > 1
            f = (q + circshift(q, (0, shift, 0)))*0.5
            i == -1 ? f[:,end,:] = q[:,end,:] : f[:,i,:] = q[:,i,:]
        elseif axis == 3 && size(q)[3] > 1
            f = (q + circshift(q, (0, 0, shift)))*0.5
            i == -1 ? f[:,:,end] = q[:,:,end] : f[:,:,i] = q[:,:,i]
        else
            f = q
        end

        return f

    end # end dnup

    function xdn(f)
        if patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger" || patch["kind"][1:7] == "ramses"
            return dnup(f, 1, 1)
        else
            return f
        end
    end # end xdn

    function ydn(f)
        if patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger" || patch["kind"][1:7] == "ramses"
            return dnup(f, 1, 2)
        else
            return f
        end
    end

    function zdn(f)
        if patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger" || patch["kind"][1:7] == "ramses"
            return dnup(f, 1, 3)
        else
            return f
        end
    end

    function xup(f)
        if patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger" || patch["kind"][1:7] == "ramses"
            return dnup(f, -1, 1)
        else
            return f
        end
    end

    function yup(f)
        if patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger" || patch["kind"][1:7] == "ramses"
            return dnup(f, -1, 2)
        else
            return f
        end
    end

    function zup(f)
        if patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger" || patch["kind"][1:7] == "ramses"
            return dnup(f, -1, 3)
        else
            return f
        end
    end

    function internal(v; all = false)

        vshape = size(v)
        if Bool(all) || length(vshape) < 3 || minimum(vshape) <= 4
            return v

        elseif patch["guard_zones"]
            l = patch["ng"]
            u = l + patch["n"]
            """ check if v.rank does not match the normal patch size.
            If so, compute the guard zone size, and adjust """

            rank = min(length(vshape), length(patch["gn"]))

            if vshape[1:rank] != tuple(patch["gn"][1:rank]...)
                gn = collect(vshape)
                ng2 = Base.copy(gn)
                for i = 1:length(patch["gn"])
                    ng2[i] = patch["gn"][i] - gn[i]
                end

                ng = ng2 .÷ 2
                n = gn .- ng2
                l = ng
                u = l + n
                #l .+= 1
            end

            return v[l[1] + 1:u[1], l[2] + 1:u[2], l[3] + 1:u[3]]

        else
            rank = min(length(vshape), length(patch["gn"]))
            if vshape[1:rank] != tuple(patch["n"][1:rank]...)
                gn = collect(vshape)
                ng2 = Base.copy(gn)

                for i = 1:length(patch["gn"])
                    ng2[i] = gn[i] - patch["n"][i]
                end

                ng = ng2 ÷ 2
                n = patch["n"]
                l = ng
                u = l + n
                return v[l[1]:u[1],l[2]:u[2],l[3]:u[3]]
            else
                return v
            end
        end
    end # end internal()

    function post_process(v; copy = false, all = false, i4 = 1)
        if copy
            v = Base.copy(v)
        end
        if ndims(v) == 4
            v = v[:, :, :, i4]
        end

        return internal(v, all = all)
    end

    """
        var(iv; all = false, copy = nothing, i4 = 1, verbose = 0)

    Evaluate arguments of various forms, including expressions.

    If the data is in spherical or cylindrical coords., then it is the angular
    momentum in the snapshot, and thus the division by metric factors.
    """
    function var(iv; all = false, copy = false, i4 = 1, verbose = 0)
        if typeof(iv) == Int
            @assert iv in patch["keys"]["numbers"] "variable index unknown"
        end

        # determine data representation
        is_ln = false
        is_vel = false
        is_hlls = false

        if strip(snap["io"]["method"]) == "parallel"
            is_ln  = ((  patch["ioformat"] - 6 )/2 == 0) ||
                     ((  patch["ioformat"] - 10)÷2 == 0)
            is_vel = ((  patch["ioformat"] - 6 )%2 == 2)

            if patch["ioformat"] == 14
                is_vel = false
                is_ln = false
            end # if
        end # if

        is_staggered = patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger"
        tot_e = patch["kind"][1:6] == "ramses"
        is_hlls = patch["kind"][1:8] == "ramses_s"

        if verbose == 1
            println("is_ln: $is_ln")
            println("is_vel: $is_vel")
            println("tot_e: $tot_e")
            println("is_staggered: $is_staggered")
        end

        v = nothing
        iv = typeof(iv) == Char ? iv = string(iv) : iv # convert to string if character

        # check if the index corresponds to a cached array
        if haskey(patch["data"], iv)
            verbose == 1 && println("$iv in data keys $(size(v))")
            v = patch["data"][iv]
        

        # or an aux item
        elseif haskey(patch["aux"]["vars"], iv)
            v = patch["aux"]["vars"][iv]["v"]
            verbose == 1 && println("$iv in keys $(size(v))")

        # or a numeric index
        elseif iv in patch["keys"]["numbers"]
            verbose == 1 && println("iv is a number: $iv")
            if iv >= 0
                v = mem(iv, verbose=verbose)
            end

        # or known variables
        elseif iv in patch["keys"]["known"]
            verbose == 1 && println("iv is a known name: $iv")

            # density
            if iv == "d"
                v = is_ln ? exp.(mem("d", verbose=verbose)) : mem("d", verbose=verbose)
            elseif iv == "lnd"
                v = is_ln ? mem("d") : log.(mem("d"))
            elseif iv == "logd"
                v = is_ln ? mem("d")/log(10) : log10.(mem("d"))

            # velocity
            elseif iv == "u1" || iv == "ux" || iv == "vx"
                if is_vel
                    v = mem("p1")
                elseif is_staggered
                    v = xup(mem("p1") ./
                        exp.(xdn(var("lnd", all=true, verbose=verbose))))
                else
                    v = mem("p1") ./ var("d", all=true, verbose=verbose)
                end 
                #v = (v)u"m/s"
            elseif iv == "u2" || iv == "uy" || iv == "vy"
                if is_vel
                    v = mem("p2")
                elseif patch["mesh_type"] != "Cartesian"
                    v = mem("p2") ./ ydn("d") ./ 
                        reshape(patch["geometric_factors"]["h2c"], :, 1, 1)
                else
                    if is_staggered
                        v = yup(mem("p2") ./
                            exp.(ydn(var("lnd", all=true, verbose=verbose))))
                    else
                        v = mem("p2") ./ var("d", all=true, verbose=verbose)
                    end
                end

            elseif iv == "u3" || iv == "uz" || iv == "vz"
                if is_vel
                    v = mem("p3")
                elseif patch["mesh_type"] != "Cartesian"
                    gf = patch["geometric_factors"]
                    v = mem("p3") ./ zdn("d") ./ 
                        reshape(gf["h31c"], :, 1, 1) ./
                        reshape(gf["h32c"], :, 1, 1)
                else
                    if is_staggered
                        v = zup(mem("p3")) ./
                            exp.(zdn(var("lnd", all=true, verbose=verbose)))
                    else
                        v = mem("p3") ./ var("d", all=true, verbose=verbose)
                    end
                end

            # Energy per unit mass
            elseif iv == "ee" || iv == "E"
                """ Expressions common to all other solvers """
                    if patch["ioformat"] == 14 || patch["ioformat"] == 15
                        v = mem("e")
                    elseif !is_staggered
                        v = (mem("e") .- 
                             var("ekin", all=true, verbose=verbose) .- 
                             var("emag", all=true, verbose=verbose)) ./ 
                             mem("d")
                    else
                        v = mem("e") ./ mem("d")
                    end

            # Kinetic energy
            elseif iv == "ekin"
                v = 0.5*var("d", all=true) .*
                    (
                        var("ux", all=true, verbose=verbose).^2 .+
                        var("uy", all=true, verbose=verbose).^2 .+
                        var("uz", all=true, verbose=verbose).^2
                    )
            # magnetic energy
            elseif iv == "emag"
                v = 0.125*(
                            xup(var("bx", all=true, verbose=verbose)) .^ 2 .+
                            yup(var("by", all=true, verbose=verbose)) .^ 2 .+
                            zup(var("bz", all=true, verbose=verbose)) .^ 2
                          )

            # Thermal energy
            elseif iv == "eth"
                if is_hlls
                    g1 = snap["gamma"] - 1.0
                    v = var("d", all=true, verbose=verbose) .^ snap["gamma"] .*
                        exp.(
                            var("s", all=true, verbose=verbose) ./ 
                            var("d", all=true, verbose=verbose)*g1)/g1
                elseif patch["gamma"] == 1.0
                    v = mem("e")
                elseif tot_e
                    v = var("e", all=true, verbose=verbose) .-
                        var("ekin", all=true, verbose=verbose) .-
                        var("emag", all=true, verbose=verbose)
                else
                    v = mem("e")
                end
            
            elseif iv == "Eth"
                v = var("eth", all=true, verbose=verbose) ./
                    var("d", all=true, verbose=verbose)

            elseif iv == "S"
                v = log.(
                         var("eth", all=true, verbose=verbose) ./
                         var("d", all=true, verbose=verbose) .^ snap["gamma"]
                         ) ./ (snap["gamma"] - 1)

            # Temperature
            elseif iv == "tt" || iv == "T"
                if is_hlls
                    g1 = snap["gamma"] - 1.0
                    vard = var("d", all=true, verbose=verbose)
                    v = vard .^ snap["gamma"] .*
                        exp.(var("s", all=true, verbose=verbose) ./
                             vard*g1) ./ vard
                else
                    v = var("Eth", all=true, verbose=verbose)
                end

            # # Magnetic field
            # elseif iv == "bx" || iv == "b1"
            #     v = xup(mem(iv))
            # elseif iv == "by" || iv == "b2"
            #     v = yup(mem(iv))
            # elseif iv == "bz" || iv == "b3"
            #     v = zup(mem(iv))
        end

        # or a letter index
        elseif iv in patch["keys"]["letters"]
            verbose == 1 && println("iv is a string: $iv $(patch["idx"]["dict"][iv])")
            iv = patch["idx"]["dict"][iv]
            if iv >= 0
                v = mem(iv)
            else
                v = 0.0*mem(1)
            end

        else
            verbose == 1 && println("unknown expression $iv, attempting to parse")
            v = evaluate_expression(patch, iv, verbose=verbose)
            return v
        end

        if v !== nothing
            """ A value v was produced, so post_process """
            return post_process(v, copy=copy, all = all, i4 = i4)
        else
            """ None of the above worked, so the iv key is not known """
            println("variable expression not understood $iv")
            return nothing
        end

    end # end var


    return var
end

""" dont know what this is """
function _h(dict)
    idx = dict["idx"]
    h = zeros(3, dict["nv"])
    if dict["kind"][1:7] == "stagger"
        if (idx["p1"] >= 1)  h[1,idx["p1"]] = -0.5 end
        if (idx["b1"] >= 1)  h[1,idx["b1"]] = -0.5 end
        if (idx["p2"] >= 1)  h[2,idx["p2"]] = -0.5 end
        if (idx["b2"] >= 1)  h[2,idx["b2"]] = -0.5 end
        if (idx["p3"] >= 1)  h[3,idx["p3"]] = -0.5 end
        if (idx["b3"] >= 1)  h[3,idx["b3"]] = -0.5 end
    end

    return h

end

""" Find the numeric index, if any, corresponding to a variable key """
function map_var(patch, iv)
    jv = nothing
    if iv in patch["keys"]["letters"]
        jv = patch["idx"]["dict"][iv]
    else
        jv = 1
    end

    return jv
end

function _pathjoin(dir, file)
    path = joinpath(dir, file)
    path = replace(path, "\\" => "/")
    return path
end