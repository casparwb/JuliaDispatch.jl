
using Printf
using JuliaDispatch.Select

""" 
    _kw_extract(kw, dict)

If key from `dict` occurs in `kw`, pop them 
"""
function _kw_extract(kw, dict)
        for (key, value) in dict
            dict[key] = pop!(kw, key, value)
        end
end

function _dir(dir, subdir)

    p = endswith(dir, '/') ? dir*subdir : dir*'/'*subdir
    if endswith(p, '/')
        return p
    else
        return p*'/'
    end
end

"""
    get_n_snapshots(;run="", data="../data")    

Returns the total number of snapshots in the given data folder
"""
function get_n_snapshots(;run="", data="../data")

    datadir = _dir(data, run)
    nsnaps = 0
    for dir in readdir(datadir)
        if isdir(datadir*dir)
            if "snapshot.nml" in readdir(datadir*dir) && !isempty(datadir*dir)
                nsnaps += 1
            end
        end
    end

    return nsnaps
end


"""
    get_snapshot_time(iout::Int; data="../data", run="")

Return the time of snapshot `iout` in the given `data/run`-folder.
"""
function get_snapshot_time(iout::Int; data="../data", run="")

    file = joinpath(_dir(data, run), "$(@sprintf("%05d", iout))")*"/snapshot.nml"
    if !isfile(file)
        println("No snapshot file with path $file")
        return nothing
    end
    time = nothing
    open(file, "r") do namelist
        for ln in eachline(namelist)
            if strip(split(ln, "=")[1]) == "TIME"
                time = parse(Float64, split(split(ln, "=")[2], ",")[1])
                break
            end
        end
    end
    return time

end

"""
    get_snapshot_ids(;data = "../data", run = "")

Return an array of IDs of all snapshots in the given `data/run`-folder.
"""
function get_snapshot_ids(;data="../data", run="")

    nsnaps = get_n_snapshots(run=run, data=data)

    datadir = _dir(data, run)
    folders = [folder for folder in readdir(datadir) if (isdir(datadir*"$folder") && 
                                                         startswith(folder, "0")  && 
                                                         "snapshot.nml" in readdir(datadir*"$folder"))]
    IDs = zeros(Int, length(folders))
    if nsnaps != length(folders) @warn "nsnaps != length(folders)" end
    for (i, folder) in enumerate(folders)
        IDs[i] = parse(Int, folder)
    end

    return IDs

end


"""
    get_new_snapshot(current_snap; data="data/", run="", sleeptime=10, maxsleep=100)

Return snapshot id with `ID = current_snap + 1`. Will wait for new snapshots as they are being produced. 
If program waits for more than `maxsleep` seconds without finding a new snapshot, the program will exit.
"""
function get_new_snapshot(current_snap; data="data/", run="", sleeptime=10, maxsleep=100)
    total_slept = 0
    while true
        snaps = get_snapshot_ids(data=data, run=run)
        new = findall(snaps .> current_snap)
        if !isempty(new)
            current_snap = snaps[new[1]]
            println("Parsing snapshot $current_snap")
            if (time() - ctime(_dir(data, run)*@sprintf("%05d", current_snap)*"/snapshot.nml")) < 3
                sleep(3)
            end             

            return current_snap#snapshot(current_snap, data=data, run=run, progress=false, suppress=true)
        else
            println("Waiting for new snapshot.")
            sleep(sleeptime)
            total_slept += sleeptime
        end

        if total_slept >= maxsleep
            println("Waited $(maxsleep)s without new snapshot. Finshing up.")
            return nothing
        end
    end
        
end


"""
    get_n_patches(iout; data="data", run="")

Return the number of patches in the given snapshot `iout` in the `data/run`-directory.
"""
function get_n_patches(iout; data="data", run="")
    snapdir = _dir(data, run)*@sprintf("%05d", iout)
    return length([dir for dir in readdir(snapdir) if endswith(dir, ".dat")])
end



function purge_cached_namelists(iout, ;data="data/", run="")
    rundir = _dir(data, run)*@sprintf("%05d/", iout)
    for entry ∈ readdir(rundir)
        if endswith(entry, ".jld2")
            rm(rundir*entry)
        end
    end
end