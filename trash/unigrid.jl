function unigrid_plane(snap::Dict; x = nothing, y = nothing, z = 0.5,
                      iv = 0, verbose=0, with_gz=false)
    """
    Function for producing a buffer of a 2d slice at the given coordinate

    Input:
    ------------------------
    - snap: dict, snapshot
    - x, y, z: float, at what position to slice
    - iv: string/int/array, what quantity(/ies) to buffer. if array, output will
        be a dictionairy with keys equal to the input variable names, and
        values as the 2d-buffers.
    - with_gz: bool, whether to include guard zones (if available)

    Output:
    -----------------------
    - dictionairy of arrays or array of floats
    """

    xyz = (x, y, z)
    ax = getindex((1, 2, 3), xyz .!= nothing)
    dims = copy(snap["cartesian"]["dims"])
    deleatat!(dims, ax)

    patches = patches_in(snap, x=x, y=y, z=z)
    if length(patches) == 0
      throw(ErrorException(" no patches found in [$x, $y, $z]"))
    end

    verbose > 0 ? println("number of patches: $(length(patches))") : nothing

    if with_gz
        ns = getindex(p0["gn"], [1, 2, 3] .!= i)
        datashp = dims .* ns
    else
        ns = getindex(p0["n"], [1, 2, 3] .!= i)
        datashp = dims .* ns
    end

    buffer = nothing
    ivs = nothing
    if typeof(iv) <: String && iv == "all"
        ivs = [iv for (k, iv) in snap["idx"] if (typeof(k) <: String && iv > 0)]
        buffer = Dict{Union{String, Int}, Array{Float32, 2}}()
        for iv in ivs
            buffer[iv] = zeros(Float32, datashp...)
        end
    elseif typeof(iv) <: Array
        ivs = iv
        buffer = Dict{Union{String, Int}, Array{Float32, 2}}()
        for iv_ in iv
            buffer[iv_] = zeros(Float32, datashp...)
        end
    else
        ivs = [iv]
        buffer = Dict{Union{String, Int}, Array{Float32, 2}}()
        buffer[iv] = zeros(Float32, datashp...)
    end

    patchDict = Dict{Int, Tuple{Array{Float32, 2}, NTuple{4, Int}, Int}}()

    for p in patches
        im = plane(p, x = x, y = y, z = z, iv = iv,
                   verbose=verbose, with_gz=with_gz)

        idxs = corner_indices(snap, p, dir=i)
        buffer[idxs[1]:idxs[2], idxs[3]:idxs[4]] = im

        patchDict[p["id"]] = (im, idxs, p["level"])

    end

    # datashp = nothing
    #
    #
    # buffer = zeros(Float32, datashp...)
    # verbose >= 1 ? println("Buffer shape: $datashp") : nothing
    #
    # sortedkeys = sort(collect(keys(ll)))
    # for id in sortedkeys
    #     im, idxs = patchDict[id][1:2]
    #     buffer[idxs[1]:idxs[2], idxs[3]:idxs[4]] = im
    # end

    return buffer
end



function unigrid_volume(snap; iv = 0, verbose=0, with_gz=false)
    """
    Function for producing a 3-dimensional volume buffer.

    Input:
    ------------------------
    - snap: dict, snapshot
    - iv: string/int/array, what quantity(/ies) to buffer. if array, output will
        be a dictionairy with keys equal to the input variable names, and
        values as the 2d-buffers.
    - with_gz: bool, whether to include guard zones (if available)

    Output:
    -----------------------
    - dictionairy of arrays or array of floats
    """

    dims = snap["cartesian"]["dims"]
    patches = snap["patches"]

    patchDict = Dict{Int, Tuple{Array{Float32, 3}, Array{Float64, 2}, Int}}()

    datasize = 0
    for p in patches
        Box = box(p, iv = iv, verbose = verbose, with_gz=with_gz)
        datasize += prod(size(Box))

        e = p["extent"]

        patchDict[p["id"]] = (Box, e, p["level"])

    end

    datashp = nothing
    try
        n = Int(cbrt(datasize))
        datashp = (n, n, n)
    catch
        if with_gz
            datashp = p0["gn"] .* dims
        else
            datashp = p0["n"] .* dims
        end
    end

    data = reshape(zeros(Float32, datasize), datashp...)
    verbose >= 1 ? println("box shape $datashp") : nothing

    sortedkeys = sort(collect(keys(ll)))
    for id in sortedkeys
        Box, e = patchDict[id][1:2]

        xs = Int.(round.(e[3,1:2] * datashp[1]))
        ys = Int.(round.(e[3,3:4] * datashp[2]))
        zs = Int.(round.(e[2,1:2] * datashp[3]))

        data[(xs[1] + 1):xs[end],
             (ys[1] + 1):ys[end],
             (zs[1] + 1):zs[end]] = Box
    end

    return data

end
