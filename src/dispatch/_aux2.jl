using FortranFiles, Printf, Mmap

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
    pos = 1
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

            tpe = String(read(file, FString{2}))

            rank = length(shp)

            verbose > 2 && println("type: $tpe")
            if startswith(tpe, "r")
                read(file, (Float32, len))
            else
                read(file, (Int32, len))
            end
        
            # v = reshape(v, Int.(shp)...)

            vars[name] = Dict("type" => tpe, "name" => name, "rank" => rank,
                                "shape" => tuple(shp...), "v" => nothing, "pos" => pos)
            pos += 1
        catch e
            break
        end

    end
    FortranFiles.close(file)
    auxDict["vars"] = vars
end

function aux_mem(name, patch)

    var = patch["aux"]["vars"][name]
    pos = var["pos"]
    shape = var["shape"]
    tpe = startswith(var["type"], "r") ? Float32 : Int32
    # data = Array{tpe, length(shape)}(undef, shape...)

    file = FortranFile(patch["aux"]["filename"])

    # offset = 2+4+4*var["rank"]+2+6
    for i = 1:pos*5
        read(file)
    end

    data = read(file, (tpe, prod(shape)))
    # data = Mmap.mmap(patch["aux"]["filename"], Array{tpe, var["rank"]}, shape, offset)
    close(file)

    data = reshape(data, Int.(shape)...)
    return data
end