using JLD

#include("../select/_select.jl")


function save_to_JLD(data::AbstractArray; var_name::String; folder=nothing,
                     filename=nothing, verbose=0)
    """
    Save the given data buffer to a JLD file, for easier
    loading later.
    """
    if var_name == nothing
        throw(ErrorException, "please give variable name")
    end

    if folder != nothing
        if !isdir(folder)
            println("creating folder $folder")
            mkdir(folder)
        end
    end

    if folder != nothing
        outpath = folder*filename*".jld"
    else
        outpath = filename*".jld"
    end

    println("saving to $outpath")
    JLD.save(outname, var_name, data)

end

function save_to_JLD(data::Dict; ivs=nothing, folder=nothing,
                     filename=nothing, verbose=0)
    """
    Save the given data buffer to a JLD file, for easier
    loading later.
    """
    if filename == nothing
        throw(ErrorException, "please give output file name")
    end

    if folder != nothing
        if !isdir(folder)
            println("creating folder $folder")
            mkdir(folder)
        end
    end


    if ivs != nothing
        if ivs <: AbstractArray
            vars = ivs
        else
            vars = [ivs]
        end
    else
        vars = collect(keys(data))
    end

    if folder != nothing
        outpath = folder*filename*".jld"
    else
        outpath = filename*".jld"
    end


    println("saving to $outpath")
    jldopen(outpath, "w") do file
        for var in vars
            var = string(var)
            JLD.write(file, var, data[var])
        end
    end


end


function load_from_JLD(var_name, folder="jld_data/"; filename=nothing)


end
