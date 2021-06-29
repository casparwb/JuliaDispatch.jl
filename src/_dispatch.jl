using PyCall, Printf, Mmap, StaticArrays

include("_aux2.jl")
include("_dispatch_grid.jl")

"""
    snapshot(iotu::Int; run::String, data::String, verbose::Int)

Parses patches and produces a `Dictionairy` with all properties of snapshot
`iout`.

Arguments:
------------
    - iout: Int, snapshot ID

Kwargs:
------------
    - data: String, path to data (snapshots), default "..\\data"
    - run: String, path to snapshots folders relative to data, default ""

Returns:
------------
    - Dictionary of snapshot
"""
function snapshot(iout=0; run="", data="..\\data", verbose = 0)

    ### wrapper for f90nml reader ###
    read_nml(file) = pyimport("f90nml").read(file)

    ### find data- and rundirs ###
    # rundir = _dir(data, run)
    statedict = Dict{String, Any}()
    rundir = joinpath(pwd(), data)
    if !isdir(rundir)
        println("directory $(rundir) does not exist")
        return nothing
    end

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

    params_list = read_nml(file)

    _add_nml_list_to(statedict, params_list)

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
                    v >= 1 ? idx["vars"][v] = k : nothing
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
    verbose > 1 ? println("checking if there is a snapshots.dat") : nothing
    try
        statedict["datafile"] = _file(rundir, "snapshots.dat")
        statedict["datafiled"] = open(statedict["datafile"], "rb")
    catch
        nothing
    end
    statedict["dict"] = statedict["nml_list"]["snapshot_nml"]

    ### dress snapshot with attributes for snapshot_nml ###
    ### lists are already arrays ###

    """ add patches as a list of dicts """
    if verbose > 1 println("add patches as a list of dicts") end

    statedict["patches"] = Array{Dict, 1}([])
    files = [f for f in readdir(datadir) if endswith(f, "_patches.nml")]

    # ignore left-over *_patches.nml from other runs
    save = []
    rank = nothing
    for i = 1:statedict["mpi_size"]
        f = files[i]
        # split name to get rank number
        rank = parse(Int, split(f, "_")[2])

        # add parameter groups from data/run/NNNNN/*patches.nml
        file = _file(datadir, f)
        if verbose > 1 println("parsing $file") end

        patch_dict = parse_patches(statedict, file)
        push!(save, patch_dict)
    end


    for patch_dict in save
        ids = sort(collect(keys(patch_dict)))
        for id in ids
            p = _patch(id, patch_dict[id], statedict, rank)
            _add_axes(statedict, p)
            push!(statedict["patches"], p)
        end
    end

    return statedict

end

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

function parse_patches(snap::Dict, file="../data/00000/rank_00000_patches.nml")

    prop_dict = Dict{Int16, Dict}()

    nml_version = 1
    nbor_info = nothing

    if "params_list" in keys(snap)
        params_list = snap[params_list]

        if "io_params" in params_list
            io_params = params_list["io_params"]

            if "nml_version" in io_params
                nml_version = io_params["nml_version"]
            end
        end
    end

    if nml_version > 1
        patches_nml = read_nml(file)
        if "nbor_nml" in keys(patches_nml)
            nblor_nml = patches_nml["nbor_nml"]
            nbor_info = Dict()
            for nb in nbor_nml
                id = nb["id"]
                parents_id = nb["parent_id"]
                nbor_ids = nb["nbor_ids"]
                nbor_info[id] = Dict("parent_id" => parent_id,
                                     "nbor_ids"  => nbor_ids )
            end
        end
    end
    open(file, "r") do fo
        watch_block = false
        d = nothing
        id = nothing
        @inbounds for line in eachline(fo)
            # strip commas and equal sign from line and split
            line = replace(
                   replace(
                   replace(line, "=" => " "), "," => ""), "''" => " ")
            items = split(line)
            if length(items) == 0
                nothing
            elseif items[1] == "&PATCH_NML"
                d = Dict{String, Any}()
                watch_block = true
            elseif items[1] == "ID"
                id = parse(Int, items[2])
            elseif items[1] == "POSITION"
                d["position"] = _parse_namelist(items)
            elseif items[1] == "SIZE"
                d["size"] = _parse_namelist(items)
                if "n" in keys(snap)
                    n = snap["n"]
                else
                    n = 1
                end
                if !("ds" in keys(d))
                    if snap["no_mans_land"]
                        d["ds"] = d["size"]./n
                    else
                        d["ds"] = d["size"]./([(snap["n"] - 1) > 1 ? (snap["n"] - 1) : 1])
                    end
                end
            elseif items[1] == "LEVEL"
                d["level"] = parse(Int32, items[2])
            elseif items[1] == "DTIME"
                d["dtime"] = parse(Float64, items[2])
            elseif items[1] == "TIME"
                d["time"] = parse(Float64, items[2])
            elseif items[1] == "ISTEP"
                d["istep"] = parse(Int32, items[2])
            elseif items[1] == "DS"
                d["ds"] = _parse_namelist(items)
            elseif items[1] == "MESH_TYPE"
                d["mesh_type"] = parse(Int32, items[2])
            elseif items[1] == "NCELL"
                d["ncell"] = _parse_namelist(items)
            elseif items[1] == "KIND"
                d["kind"] = items[2]
            elseif items[1] == "/" && watch_block # the final entry of a namelist is always "/"
                if nbor_info != nothing
                    d["parent_id"]=nbor_info[id]["parent_id"]
                    d["nbor_ids"]=nbor_info[id]["nbor_ids"]
                end
                prop_dict[id]=d
                watch_block = false
            end
        end # end eachline
    end # end open
    return prop_dict
end # end parse_patches

function _parse_namelist(items)

    pos = items[2:end]
    if length(pos) == 1
        pos[1] = replace(pos[1], "3*" => "")
        pos = [pos[1], pos[1], pos[1]]
    elseif length(pos) == 2
        if findfirst("2*", pos[1]) != nothing
            pos[1] = replace(pos[1], "2*" => "")
            pos = [pos[1], pos[1], pos[2]]
        end
        if findfirst("2*", pos[2]) != nothing
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

    #endswith(dir, "/") ? p = dir*subdir : dir*"/"*subdir
    if endswith(dir, "\\")
        p = dir*subdir
    else
        p = dir*"\\"*subdir
    end
    endswith(p, "\\") ? p : p*"\\"
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


function _patch(id, patch_dict, snap, rank, verbose=0)

    patch = Dict()
    patch["id"] = id
    patch["rank"] = rank
    patch["memmap"] = 1

    # add general attributes from snapshot.nml
    for (k, v) in snap["dict"]
        patch[k] = v
    end

    for (k, v) in patch_dict
        patch[k] = v
    end

    if !patch["guard_zones"]
        patch["li"][:] .= 1
        patch["ui"][:] .= patch["n"]
    end

    if "idx" in keys(snap)
        patch["idx"] = snap["idx"]
        patch["idx"]["h"] = _h(patch)
    end

    # reconstruct items pruned from patch_nml

    if "size" in keys(patch) && "position" in keys(patch)
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

    if !("units" in keys(patch))
        if "units" in keys(snap)
            patch["units"] = snap["units"]
            if !("u" in keys(patch["units"]))
                patch["units"]["u"] = patch["units"]["l"]/patch["units"]["t"]
            end
            if !("d" in keys(patch["units"]))
                patch["units"]["d"] = patch["units"]["m"]/patch["units"]["l"]^3
            end
            if !("p" in keys(patch["units"]))
                patch["units"]["p"] = patch["units"]["d"]/patch["units"]["u"]^2
            end
            if !("e" in keys(patch["units"]))
                patch["units"]["e"] = patch["units"]["m"]/patch["units"]["u"]^2
            end
            if !("b" in keys(patch["units"]))
                patch["units"]["b"] = patch["units"]["u"]*sqrt(4π*patch["units"]["d"])
            end
        end

    end

    # modify `mesh_type` from integer to string for readability
    if "mesh_type" in keys(patch)
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
                        "/$(@sprintf("%05d/%05d", patch["iout"], patch["id"])).dat"
        patch["var"] = _var(patch, patch["filename"], snap)
    else
        patch["var"] = _var(patch, snap["datafiled"], snap)
    end

    """ add a comprehensive set of variable keys """
    patch["aux"] = Dict()
    patch["data"] = Dict()
    patch["keys"] = Dict()
    patch["keys"]["letters"] = collect(keys(snap["idx"]))
    patch["keys"]["numbers"] = collect(values(snap["idx"]))
    patch["keys"]["known"] = ["d","lnd","logd","ux","uy","uz","u1","u2","u3",
                             "ee","E","T","eth","ekin"]


    # attach an aux filename, if one exists for this task
    auxfile = "/$(@sprintf("%05d", id)).aux"

    auxfile = snap["datadir"]*auxfile
    patch["auxfile"] = auxfile
    patch["aux"] = aux(id = id, rundir=snap["rundir"], datadir=snap["datadir"],
                       io = snap["iout"], file = auxfile, verbose=verbose)
    patch["keys"]["aux"] = collect(keys(patch["aux"]["vars"]))

    """ Collect all keys in a single array """
    all = []
    for key_list in values(patch["keys"])
        append!(all, key_list)
    end
    patch["all_keys"] = all

    return patch
end


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
        offset *= bytes
        push!(patch["offset"], offset)
    end

    function mem(iv; verbose = 0)
        """
        Translate alphabetic variable keys to numeric
        """

        verbose == 1 ? println("mem($iv)") : nothing

        if typeof(iv) == typeof("d")
            iv = patch["idx"][iv]
        end

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
        if patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger"
            return dnup(f, 1, 0)
        else
            return f
        end
    end # end xdn

    function ydn(f)
        if patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger"
            return dnup(f, 1, 1)
        else
            return f
        end
    end

    function zdn(f)
        if patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger"
            return dnup(f, 1, 2)
        else
            return f
        end
    end

    function xup(f)
        if patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger"
            return dnup(f, -1 ,0)
        else
            return f
        end
    end

    function yup(f)
        if patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger"
            return dnup(f, -1 ,1)
        else
            return f
        end
    end

    function zup(f)
        if patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger"
            return dnup(f, -1 ,2)
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

            rank = min(length(vshape), length(patch["ng"]))

            if vshape[1:rank] != tuple(patch["gn"][1:rank])
                gn = collect(vshape)
                ng2 = Base.copy(gn)
                for i = 1:length(patch["gn"])
                    ng2[i] = patch["gn"][i] - gn[i]
                end

                ng = ng2 .÷ 2
                n = gn .- ng2
                l = ng
                u = l + n
                l .+= 1
            end

            return v[l[1]:u[1], l[2]:u[2], l[3]:u[3]]

        else
            rank = min(length(vshape), length(patch["gn"]))
            if vshape[1:rank] != tuple(patch["n"][1:rank])
                gn = Array(vshape)
                ng2 = copy(gn)

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
        if ndims(v) == 4
            v = v[:, :, :, i4]
        end

        return internal(v, all = all)
    end

    function var(iv; all = false, copy = nothing, i4 = 1, verbose = 0)
        """
        Evaluate arguments of varuious forms, including expressions.

        If the data is in spherical or cylindrical coords., then it is the angular
        momentum in the snapshot, and thus the division by metric factors.
        """
        if typeof(iv) == Int
            @assert iv in patch["keys"]["numbers"] "variable index unknown"
        end


        # determine data representation

        is_ln = ((  patch["ioformat"] - 6 )/2 == 0
                || (patch["ioformat"] - 10)÷2 == 0
                ||  patch["ioformat"] == 15
                ||  patch["ioformat"] == 14)

        is_v  = (( patch["ioformat"] - 6)%2==0
                || patch["ioformat"] == 15
                || patch["ioformat"] == 14)

        is_st = patch["kind"][1:4] == "zeus" || patch["kind"][1:7] == "stagger"
        v = nothing

        # check if the index corresponds to a cached array
        if iv in keys(patch["data"])
            v = patch["data"][iv]

        # or an aux item
        elseif iv in keys(patch["aux"]["vars"])
            v = patch["aux"]["vars"][iv]["v"]
            verbose == 1 ? println("$iv in keys $(size(v))") : nothing

        # or a numeric index
        elseif iv in patch["keys"]["numbers"]
            verbose == 1 ? println("iv is a number: $iv") : nothing
            if iv >= 1
                v = mem(iv, verbose=verbose)
            end

        # or known variables
        elseif iv in patch["keys"]["known"]
            verbose == 1 ? println("iv is a known name $iv") : nothing

            # density
            if iv == "d"
                v = is_ln ? exp.(mem("d")) : mem("d")
            elseif iv == "lnd"
                v = is_ln ? mem("d") : log.(mem("d"))
            elseif iv == "logd"
                v = is_ln ? mem("d")/log(10) : log10.(mem("d"))

            # velocity
        elseif iv == "u1" || iv == "ux" iv == "vx"
                if is_v
                    v = mem("p1")
                elseif is_st
                    v = xup(mem("p1")) ./
                        exp.(xdn(var("lnd", all=true, verbose=verbose)))
                else
                    v = mem("p1")/var("d", all=true, verbose=verbose)
                end

            elseif iv == "u2" || iv == "uy" iv == "vy"
                if is_v
                    v = mem("p2")
                elseif patch["mesh_type"] != "Cartesian"
                    shp = size(patch["geometric_factors"]["h2c"])
                    v = mem("p2") / ydn("d") ./
                        reshape(patch["geometric_factors"]["h2c"], (shp[1], 1, 1))
                else
                    if is_st
                        v = yup(mem("p2")) ./
                            exp.(ydn(var("lnd", all=true, verbose=verbose)))
                    else
                        v = mem("p2") ./ var("d", all=true, verbose=verbose)
                    end

                end

            elseif iv == "u3" || iv == "uz" iv == "vz"
                if is_v
                    v = mem("p3")
                elseif patch["mesh_type"] != "Cartesian"
                    gf = patch["geometric_factors"]
                    shp1 = size(gf["h31c"])
                    shp2 = size(gf["h32c"])
                    v = mem("p3") / zdn("d") / reshape(gf["h31c"], (shp1[1], 1, 1)) /
                        reshape(gf["h32c"], (1, shp2[2], 1))

                else
                    if is_st
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

            # Thermal energy
            elseif iv == "eth"
                if patch["ioformat"] == 14 || patch["ioformat"] == 15 || is_st
                    v = mem("e")
                else
                    v = mem("e") .- var("ekin", all = true, verbose = verbose)
                end

            # Temperature
            elseif iv == "tt" || "T"
                v = var("eth", all = true, verbose = verbose)
                if "scaling" in keys(snap)
                    v = v*snap["scaling"]["temp"]
                end

            # Magnetic field
            elseif iv == "bx" || iv == "b1"
                v = xup(mem(iv))
            elseif iv == "by" || iv == "b2"
                v = yup(mem(iv))
            elseif iv == "bz" || iv == "b3"
                v = zup(mem(iv))

        end

        elseif iv in patch["keys"]["letters"]
            verbose == 1 ? println("iv is a letter") : nothing
            iv = patch["idx"][iv]
            v = mem(iv)

        else
            verbose == 1 ? println("unknown expression $iv") : nothing
            v = expression_parser(patch, iv) # NOT FINISHED TO-DO
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
        idx["p1"] >= 1 ? h[1,idx["p1"]] = -0.5 : nothing
        idx["b1"] >= 1 ? h[1,idx["b1"]] = -0.5 : nothing
        idx["p2"] >= 1 ? h[2,idx["p2"]] = -0.5 : nothing
        idx["b2"] >= 1 ? h[2,idx["b2"]] = -0.5 : nothing
        idx["p3"] >= 1 ? h[3,idx["p3"]] = -0.5 : nothing
        idx["b3"] >= 1 ? h[3,idx["b3"]] = -0.5 : nothing
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

""" Returns which patches are in a given plane """
function patches_in(snap; x = nothing, y = nothing, z = nothing)
    patches = snap["patches"]

    if x != nothing
        patches = [p for p in patches
                   if (x >= p["extent"][3, 1] && x < p["extent"][3, 2])]
    end

    if y != nothing
        patches = [p for p in patches
                   if (y >= p["extent"][1, 1] && y < p["extent"][1, 2])]
    end

    if z != nothing
        patches = [p for p in patches
                   if (z >= p["extent"][2, 1] && z < p["extent"][2, 2])]
    end

    return patches
end

"""
plane(patch::Dict; x::Float, y::Float, z::Float,
      iv::Union{String, Int}, all::Bool, verbose::Int)

Return patch data of quantity iv at a slice x/y/z.

Arguments:
--------------
    - patch: Dictionairy, a patch object from a snapshot

Kwargs:
-------------
    - x, y, z: Float, position at which to slice, default nothing
    - iv: String/Int, what quantity to extract, default 0
    - all: Bool, whether to include guard zones, default false

Returns:
-------------
    - 2d array of Float32.
"""
function plane(patch; x = nothing, y = nothing, z = nothing, iv = 0,
                       verbose = 0, all=false)



    if patch["guard_zones"]# && !all
        li = patch["li"]  # lower inner
        ui = patch["ui"]  # upper inner
    else #all || !patch["guard_zones"]
        li = ones(Int, 3) # using Static??
        ui = patch["gn"]
    end

    if x !== nothing
        p = (x - patch["x"][1])/patch["ds"][1]
        i = Int(round(p))
        i = min(i, ui[1]-1)
        p -= i

        f = patch["var"](iv)[i, li[2]:ui[2], li[3]:ui[3]]*(1.0 - p) +
            patch["var"](iv)[i+1, li[2]:ui[2], li[3]:ui[3]]*p

    elseif y !== nothing
        p = (y - patch["y"][1])/patch["ds"][2]
        i = Int(round(p))
        i = min(i, ui[2]-1)
        p = p - i
        f = transpose(patch["var"](iv)[li[1]:ui[1], i  , li[3]:ui[3]]*(1.0-p) +
                      patch["var"](iv)[li[1]:ui[1], i+1, li[3]:ui[3]]*p)

    elseif z !== nothing
        p = (z - patch["z"][1])/patch["ds"][3]
        i = Int(round(p))
        i = min(i, ui[3]-1)
        p = p - i
        if i == 0
            f = patch["var"](iv)[li[1]:ui[1], li[2]:ui[2], end]*(1.0 - p) +
                patch["var"](iv)[li[1]:ui[1], li[2]:ui[2], i+1]*p
        else
            f = patch["var"](iv)[li[1]:ui[1], li[2]:ui[2], i]*(1.0 - p) +
                patch["var"](iv)[li[1]:ui[1], li[2]:ui[2], i+1]*p
        end
    else
        throw(ArgumentError("plane: must give one x, y, or z-value."))
    end

    Bool(verbose) ? println("plane: $i, p = $i $p") : nothing

    return f

end


"""
box(patch::Dict; iv::Union{String, Int}, all::Bool, verbose::Int)

Return volume of data of quantity iv from patch.

Arguments:
------------
    - patch: Dictionairy, a patch object from a snapshot

Kwargs:
--------------
    - iv: String/Int, what quantity to extract, default 0
    - all: Bool, whether to include guard zones, default false

Returns:
------------
    - 3d array of Float32.
"""
function box(patch; iv = 0, all=false, verbose = 0)


    if patch["guard_zones"] && !all
        li = patch["li"]  # lower inner
        ui = patch["ui"]  # upper inner
    elseif all || !patch["guard_zones"]
        li = ones(Int, 3)
        ui = patch["gn"]
    end

    f = patch["var"](iv)[li[1]:ui[1], li[2]:ui[2], li[3]:ui[3]]

    return f
end
