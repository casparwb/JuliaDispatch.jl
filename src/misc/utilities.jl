

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
