
using Printf

function _kw_extract(kw, dict)
    """ if key from dict occur in kw, pop them """
        for (key, value) in dict
            dict[key] = pop!(kw, key, value)
        end
    # return kw, dict
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
Returns the total number of snapshots in the given data folder
"""
function get_n_snapshots(;run="", data="../data")

    datadir = _dir(data, run)
    nsnaps = 0
    for dir in readdir(datadir)
        if isdir(datadir*dir)
            if "snapshot.nml" in readdir(datadir*dir)
                nsnaps += 1
            end
        end
    end

    return nsnaps
end



function get_available_ivs(snap)
    """ Function for printing what quantities
    are available in the snapshot """

    quantDict = Dict{String, String}(
                "d"    => "density",
                "lnd"  => "log density",
                "logd" => "log10 density",
                "u1"   => "x-velocity",
                "ux"   => "x-velocity",
                "vx"   => "x-velocity",
                "u2"   => "y-velocity",
                "uy"   => "y-velocity",
                "vy"   => "y-velocity",
                "u3"   => "z-velocity",
                "uz"   => "z-velocity",
                "vz"   => "z-velocity",
                "ee"   => "specific",
                "E"    => "specific",
                "ekin" => "kinetic energy",
                "eth"  => "thermal energy",
                "tt"   => "temperature",
                "T"    => "temperature",
                "bx"   => "x-magnetic field",
                "by"   => "y-magnetic field",
                "bz"   => "z-magnetic field",
                "b1"   => "x-magnetic field",
                "b2"   => "y-magnetic field",
                "b3"   => "z-magnetic field",
                "px"   => "x-momentum",
                "py"   => "y-momentum",
                "pz"   => "z-momentum",
                "p1"   => "x-momentum",
                "p2"   => "y-momentum",
                "p3"   => "z-momentum",
                "e"    => "energy",
                "s"    => "",
                "ex"   => "x-electric field",
                "ey"   => "y-electric field",
                "ez"   => "z-electric field"
                )


    for (key, val) in snap["idx"]
        if typeof(val) <: Int && val > 0
            try
                println("$key: $(quantDict[key])")
            catch
                println("key $key not found")
            end
        end
    end
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
    IDs = zeros(Int, nsnaps)

    datadir = _dir(data, run)
    folders = [folder for folder in readdir(datadir) if (isdir(datadir*"$folder") && 
                                                         startswith(folder, "0")  && 
                                                         "snapshot.nml" in readdir(datadir*"$folder"))]
    
    for (i, folder) in enumerate(folders)
        IDs[i] = parse(Int, folder)
    end

    return IDs

end

