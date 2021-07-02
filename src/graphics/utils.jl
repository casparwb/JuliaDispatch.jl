module Utils

export _kw_extract, get_n_snapshots, get_unit, get_available_ivs

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
    #     if "snapshot.nml" in readdir(dir)
    #         nsnaps += 1
    #     end
    # end

    # nsnaps = length([dir for dir in readdir(datadir) if startswith(dir, "0")])
    nsnaps
end



function get_unit(snap, iv)#, system=nothing)

    # if system == nothing
    if !haskey(snap, "units") return nothing end
    system = strip(snap["units"]["system"])
    # end

    system = (1, 2)[system .== ["si", "cgs"]][1]
    # index 1 is SI and index 2 is CGS

    units = Dict{String, Array{String, 1}}(
            "dens"          => [" kg/m^3", "g/cm^3"],
            "vel"        => ["m/s", "cm/s"],
            "energypm" => ["J/kg", "\\frac{g cm^2}{s^2 g}"],
            "energy"          => ["J", "g cm^2 s^{-2}"],
            "temp"       => ["K"],
            "mag"    => ["T", "G"],
            "momentum"   => ["kg m s^{-1}", "g cm s^{-1}"]
            )

    unit = nothing
    quant = nothing
    if any(iv .== ("d", "lnd", "logd"))
        unit = units["dens"][system]
        quant = "density "
    elseif any(iv .== ("u1", "u2", "u3", "ux", "uy", "uz", "vx", "vy", "vz"))
        unit = units["vel"][system]
        quant = "velocity "
    elseif any(iv .== ("tt", "T"))
        unit = units["temp"][1]
        quant = "temperature "
    elseif iv == "ekin"
        unit = units["energy"][system]
        quant = "kinetic energy "
    elseif iv == "eth"
        unit = units["energy"][system]
        quant = "thermal energy "
    elseif any(iv .== ("ee", "E"))
        unit = units["energypm"][system]
        quant = "specific energy "
    elseif any(iv .== ("b1", "b2", "b3", "bx", "by", "bz"))
        unit = units["mag"][system]
        quant = "\\text{magnetic field} "
    end

    unit = "\$ \\left[" *unit* "\\right] \$"

    # return latexstring(quant*unit)
    return latexstring(quant*unit)

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

end #module