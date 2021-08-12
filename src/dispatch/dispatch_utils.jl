
using Printf
using JLD2
using JuliaDispatch.Utils



"""
    read_patch_metadata(io = 0, run = "", data="../data", n_rank = 0, verbose = 0)

Read patch metadata from Fortran unformatted file
"""
function read_patch_metadata(io = 0, run = "", data="../data", n_rank = 0; verbose = 0)
    rundir = _dir(data, run)
    if iszero(n_rank)
        dir = rundir*@sprintf("%05d", io)
        files = [file for file in readdir(dir) if endswith(file, "patches.nml")]
        n_rank = length(files)
    end # if

    patches = Dict()
    file0 = rundir*@sprintf("%05d", io)*"/patches.jld2"
    file1 = rundir*@sprintf("%05d", io)*"/rank_00000_patches.nml"
    if ispath(file0) && ctime(file0) > ctime(file1)
        patches = JLD2.load(file0)
        verbose > 0 && @info "$file0, $(length(keys(patches))) patches"
    else
        verbose > 0 && @info "Caching metadata to $file0"

        rank = 0
        while rank < n_rank
            file1 = rundir*@sprintf("%05d", io)*"/rank_"*@sprintf("%05d", rank)*"_patches.nml"
            read_rank_namelist(file1, patches, rank, verbose=verbose) 
            verbose > 1 && println("$file1, $(length(keys(patches))) patches")
            rank += 1
        end # while

        try
            JLD2.save(file0, patches)
        catch e
            throw(e)
        end # try
    end # if

    np = length(keys(patches))
    verbose > 0 && @info "$n_rank rank(s), $np patches"

    if iszero(np)
        println("no files found")
        return nothing
    else 
        return patches
    end # if  
end # function

"""
    read_rank_namelist(file, patches, rank; verbose=0)

Optimized parsing of a patch namelist entries. Updates dictionary with
patch ID as key and a dictionary with multi properties as the value
"""
function read_rank_namelist(file, patches, rank; verbose=0)
    n = 0
    
    open(file, "r") do fo
        nbor_block = false
        watch_block = false
        
        # the following two assinments are require for them to remain in scope
        d = nothing
        id = nothing
        nbor_ids = nothing
        for line in eachline(fo)
            # strip commas and equal sign from line and split
            line = replace(replace(replace(line, "=" => " "), 
                                                    "," => ""), 
                                                    "''" => " ")
            items = split(line)
            if length(items) == 0
                continue
            end # if
            #ds = []

            # detect blocks to parse
            if items[1] == "&PATCH_NML"
                d = Dict()
                watch_block = true
            elseif items[1] == "&NBOR_NML"
                nbor_block = true
                verbose > 0 && println("starting nbor block")
            
            # parse a block of the namelist with an nbor list
            elseif nbor_block
                function parse_nbors!(items, nbor_ids)
                    for item in items
                        if occursin('*', item)
                            println("warning; multiple $id $item")
                            rep = split(item, '*')
                            for i = 1:parse(Int, rep[1])
                                push!(nbor_ids, parse(int, rep[1]))
                            end # for
                        else
                            push!(nbor_ids, parse(Int, item))
                        end # if
                    end # for
                end # function

                verbose > 1 && println(items)

                itm0 = items[1]

                if itm0 == "ID"
                    id = strip(items[2])#parse(Int, items[2])
                    verbose > 0 && println("nbor_block id: $id")
                elseif itm0 == "PARENT_ID"
                    parent_id = parse(Int, items[2])
                    verbose > 0 && println("nbor_block parent_id: $parent_id")
                elseif itm0 == "NBOR_IDS"
                    nbor_ids = []
                    parse_nbors!(items[2:end], nbor_ids)
                    verbose > 0 && println("nbor_block nbor_ids: $nbor_ids")
                elseif string(itm0) == "/"
                    nbor_ids = unique(sort(nbor_ids))
                    patches[id]["nbor_ids"] = nbor_ids
                    nbor_block = false
                    verbose > 0 && println("ending nbor_block")
                else
                    parse_nbors!(items, nbor_ids)
                    verbose > 0 && println("appended nbor_ids: $items")
                end # if

            elseif watch_block
                length(items) == 0 && continue
                itm0 = items[1]
                if itm0 == "ID"
                    id = parse(Int, items[2])
                elseif itm0 == "POSITION"
                    d["position"] = _parse_multiple(items)
                elseif itm0 == "SIZE"
                    d["size"] = _parse_multiple(items)
                elseif itm0 == "LEVEL"
                    d["level"] = parse(Int, items[2])
                elseif itm0 == "DTIME"
                    d["dtime"] = parse(Float64, items[2])
                elseif itm0 == "TIME"
                    d["time"] = parse(Float64, items[2])
                elseif itm0 == "ISTEP"
                    d["istep"] = parse(Int, items[2])
                elseif itm0 == "RECORD"
                    d["record"] = parse(Int, items[2])     
                elseif itm0 == "DS"
                    d["ds"] = _parse_multiple(items)
                elseif itm0 == "MESH_TYPE"
                    d["mesh_type"] = parse(Int, items[2]) 
                elseif itm0 == "NCELL"
                    d["ncell"] = _parse_multiple(items)
                elseif itm0 == "KIND"
                    d["kind"] = items[2]
                elseif string(itm0) == "/" && watch_block # the final entry is always "/"
                    d["rank"] = rank
                    patches[string(id)] = d
                    watch_block = false
                    n += 1
                    verbose > 2 && println("id $id")           
                end # if
            
            end # if

        end # for
            

    end # do

    verbose > 1 && println("$file: $n patches")

end # function                

"""
    _parse_multiple(items)

Figure out namelist records that may have N*xxx syntax
"""
function _parse_multiple(items)
    pos = items[2:end]
    if length(pos) == 1
        pos[1] = replace(pos[1], "3*" => "")
        pos = [pos[1], pos[1], pos[1]]
    elseif length(pos) == 2
        if occursin("2*", pos[1])
            pos[1] = replace(pos[1], "2*" => "")
            pos = [pos[1], pos[1], pos[2]]
        end # if
        if occursin("2*", pos[2])
            pos[2] = replace(pos[2], "2*" => "")
            pos = [pos[1], pos[2], pos[2]]
        end #if
    end # if
    try 
        pos = [parse(Int, p) for p in pos]
    catch
        pos = [parse(Float64, p) for p in pos]
    end # try

    return pos
end # function


"""
    cache_snapshots_live(;data="data/", run="", current_snap=0, sleeptime=10)

Read snapshots live as they are being produced by a simulation in the given `data/run` directory, caching the 
namelists for faster loading of snapshots later.
 Will start from snapshot `current_snap`.
 Will wait for `sleeptime` seconds to look again if no new snapshot is found. If program waits for more than
`maxsleep` seconds without finding a new snapshot, the program will exit.
"""
function cache_snapshots_live(;data="data/", run="", current_snap=0, sleeptime=10, maxsleep=100)
    while true
        snap = get_new_snapshot(current_snap, data=data, run=run, sleeptime=sleeptime, 
                                              maxsleep=maxsleep)
        if isnothing(snap)
            break
        else
            current_snap = snapshot(snap, data=data, run=run, progress=false)["iout"]
        end
    end
    return nothing
end



function get_all_snapshots(;data=data, run=run)

    snapIDs = get_snapshot_ids(data=data, run=run)
    snapshots = Vector{Dict}(undef, length(snapIDs))

    for idx in ProgressBar(eachindex(snapIDs))
        snapshots[idx] = snapshot(snapIDs[idx], data=data, run=run, progress=false, verbose=-1)  
    end

    return snapshots
end


"""
    get_snapshots(;data=data, run=run, tspan=nothing)

Returns a vector of parsed snapshots in the given  `data/run` folder. If `tspan` is given, only
snapshots with `tspan[1] <= "time" <= tspan[2]` will be parsed and returned.
"""
function get_snapshots(;data=data, run=run, tspan=nothing)

    if isnothing(tspan) 
        return get_all_snapshots(data=data, run=run) 
    else
        start, stop = tspan
        snapIDs = get_snapshot_ids(data=data, run=run)
        times = get_snapshot_time(snapIDs, data=data, run=run)

        indices = findall(start .<= times .<= stop)

        snapIDs = snapIDs[indices]

        snapshots = Vector{Dict}(undef, length(snapIDs))

        for idx in ProgressBar(eachindex(snapIDs))
            snapshots[idx] = snapshot(snapIDs[idx], data=data, run=run, progress=false, verbose=-1)  
        end

        return snapshots
    end
end