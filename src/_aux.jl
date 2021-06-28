import Conda

#ffreader = pyimport("scipy.io").FortranFile

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
        aux_read(auxDict, verbose=verbose)
    end
    return auxDict
end

function aux_read(auxDict; verbose=0)
    ffreader = pyimport_conda("scipy.io", "scipy").FortranFile
    ff = ffreader(auxDict["filename"], "r")
    auxDict["version"], auxDict["id"] = ff.read_ints()

    vars = Dict()
    while true
        try
            ii = ff.read_record("b")
            name = to_str(ii)

            verbose > 2 ? println(" name: $name") : nothing

            rnk = ff.read_ints()[1]

            verbose > 2 ? println(" rank: $rnk") : nothing

            shp = ff.read_ints()

            verbose > 2 ? println(" shape: $shp") : nothing

            rank = length(shp)

            ii = ff.read_record("b")
            tpe = to_str(ii)

            verbose > 2 ? println("type: $tpe") : nothing
            if string(tpe[1]) == "r"
                v = ff.read_reals("<f4")
            else
                v = ff.read_ints("<i4")
            end
            v = reshape(v, reverse(shp)...)
            if ndims(v) == 2
                v = transpose(v)
            else
                # v = permutedims(v, (3, 2, 1))
            end
            vars[name] = Dict("type" => tpe, "name" => name, "rank" => rank,
                                "shape" => shp, "v" => v)
        catch
            break
        end

    end
    ff.close()
    auxDict["vars"] = vars
end

function to_str(ii)
    s = ""
    for i in ii
        i != 32 ? s *= string(Char(i)) : nothing
    end

    return s


end
