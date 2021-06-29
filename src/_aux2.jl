using FortranFiles

function aux(;id = 1, io = 0, rundir = "", datadir = "../data", file=nothing,
    verbose = 0)
    """ aux format file object """
    auxDict = Dict()
    auxDict["vars"] = Dict()

    if isfile(file)
        auxDict["filename"] = file
    else
        auxDict["filename"] = datadir*"/"*rundir*"/$(@sprintf("%05d/%05d", io, id)).aux"
    end
    if isfile(file)
        aux_read!(auxDict, verbose=verbose)
    end
    return auxDict
end

function aux_read!(auxDict; verbose=0)

    file = FortranFile(auxDict["filename"])
    auxDict["version"], auxDict["id"] = read(file, (Int32, 2))

    vars = Dict()
    while true
        try
            name = String(read(file, FString{2}))

            verbose > 2 && println(" name: $name")

            rnk = read(file, Int32)

            verbose > 2 && println(" rank: $rnk")
            shp = read(file, (Int32, rnk))

            verbose > 2 && println(" shape: $shp")

            rank = length(shp)
            len = prod(shp)
            tpe = read(file, FString{2})

            rank = length(shp)

            verbose > 2 && println("type: $tpe")

            if startswith(String(tpe), "r")
                v = read(file, (Float32, len))
            else
                v = read(file, (Int32, len))
            end
        
            v = reshape(v, Int.(shp)...)

            vars[name] = Dict("type" => tpe, "name" => name, "rank" => rank,
                                "shape" => shp, "v" => v)
        catch e
            break
        end

    end
    FortranFiles.close(file)
    auxDict["vars"] = vars
end
