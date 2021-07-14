

using JuliaDispatch.Select, Unitful
import Interpolations
const Itp = Interpolations
"""
    init_buffer(snap::Dict, iv::Union{Int, String}, dims::Union{Tuple, AbstractArray},
                ndims::Int)

Initialize a buffer for storing the data. Returns a `Dict` of Arraya if `iv` is a collection, otherwise
an `Array`.

#Arguments:
- `snap::Dict`: snapshot

#Kwargs:
- `iv::Union{String, Int, Collection}` of String/Ints, what quantity(/ies) to extract, default 0
- `dims::Union{Int, Tuple, Array}`, data size in each dimension. If Int,
                                    all dimensions will have same length. If Tuple/Array,
                                    must have length(dims) == ndims.
- `ndims::Int`: number of dimensions
"""
function init_buffer(snap, iv, dims, num_dims)

    datashp = nothing
    # determine datashape 
    if typeof(dims) <: Int
        datashp = repeat([dims], num_dims) # same number of points in each dimension
    else
        datashp = dims # or array of points in each dimension
    end

    buffer = Dict{Union{String, Int}, Array{Number, num_dims}}() # dict for storing quantities

    if typeof(iv) <: String && iv == "all" # if the given iv is `all`
        ivs = [k for (k, iv) in keys(snap["idx"]["dict"])
               if (typeof(iv) <: Int && iv > 0)]

        for iv in ivs
            # buffer[iv] = Array{Number, num_dims}(0, datashp...)
            buffer[iv] = zeros(Float32, datashp...)
        end
    elseif typeof(iv) <: AbstractArray # if the given iv is an array of different ivs
        for iv_ in iv
            # buffer[iv_] = Array{Number, num_dims}(0, datashp...)
            buffer[iv_] = zeros(Float32, datashp...)
        end
    else # if the iv is a single entry
        # buffer[iv] = Array{Number, num_dims}(0, datashp...)
        buffer[iv] = zeros(Float32, datashp...)
    end

    return buffer

end

"""
    amr_plane(snap::Dict; iv::Union{Int, String}, x::Float, y::Float, z::Float,
             Log::Bool, dims::Union{Int, Tuple})

Return a 2D array containing interpolated data of quantity `iv` in a slice `x/y/z` from all patches in a given snapshot. If `dims` is an
`Int` the resulting array will be of size `(dims, dims)`. If `dims` is a length-2 array, the array will have size `(dims...)`. 
"""
function amr_plane(snap; iv = 0, x = nothing, y = nothing, z = nothing,
                   Log = false, dims::Union{Int, Tuple}=100, all=false)


    xyz = [x, y, z]
    axIdx = getindex((1, 2, 3), xyz .!= nothing)[1]
    ax = xyz[axIdx]
    dir1, dir2 = getindex((1, 2, 3), xyz .== nothing)
    dir1s, dir2s = getindex(("x", "y", "z"), xyz .== nothing)


    # e1s, e2s = nothing, nothing
    # if axIdx == 1 || axIdx == 2
    #     e1s = [1, 2]
    #     e2s = [3, 4]
    # elseif axIdx == 3
    #     e1s = [3, 4]
    #     e2s = [1, 2]
    # end

    Size = copy(snap["cartesian"]["size"])
    origin = copy(snap["cartesian"]["origin"])

    deleteat!(Size, axIdx)
    deleteat!(origin, axIdx)

    buffer = init_buffer(snap, iv, dims, 2)
    if typeof(dims) <: Int
        n1, n2 = dims, dims
    else
        n1, n2  = dims
    end

    patches = patches_in(snap, x=x, y=y, z=z)
    if length(patches) == 0
      throw(ErrorException(" no patches found in [$x, $y, $z]"))
    end

    d1 = range(origin[1], origin[1]+Size[1], length=n1) # axis in plane dimension 1
    d2 = range(origin[2], origin[2]+Size[2], length=n2) # axis in plane dimension 2

    interpx(itp, y, z) = itp(ax, y, z)
    interpy(itp, x, z) = itp(x, ax, z)
    interpz(itp, x, y) = itp(x, y, ax)
    interp = [interpx, interpy, interpz][axIdx] 

    for iv in keys(buffer)
        for patchID in eachindex(patches)#patch in patches
            patch = patches[patchID]
            # interpolate between planes
            # data = plane(patch, iv=iv, x=x, y=y, z=z, all=true)


            d1extent = patch[dir1s]
            d2extent = patch[dir2s]

            d1s = findall(d1extent[1] .< d1 .< d1extent[end]) # indices extended by patch in plane dimension 1
            d2s = findall(d2extent[1] .< d2 .< d2extent[end]) # indices extended by patch in plane dimension 2

            # interpolate patch data
            # itp = Itp.interpolate((patch[dir1s], patch[dir2s]), data, Itp.Gridded(Itp.Linear()))
            itp = Itp.interpolate((patch["x"], patch["y"], patch["z"]), box(patch, iv=iv, all=true), Itp.Gridded(Itp.Linear()))
            for i2 ∈ d2s, i1 ∈ d1s
                # buffer[iv][i1, i2] = itp(d1[i1], d2[i2])
                buffer[iv][i1, i2] = interp(itp, d1[i1], d2[i2]) # insert into buffer
            end
        end
    end

    if length(keys(buffer)) == 1
        return collect(values(buffer))[1]
    else
        return buffer
    end

end


"""
    amr_volume(snap::Dict; iv::Union{Int, AbstractArray, String}, all::Bool,
            dims::Union{Tuple, Int}, verbose::int)

Return a 3D array containing interpolated data of quantity `iv` from all patches in a given snapshot. If `dims` is an
`Int` the resulting array will be of size `(dims, dims, dims)`. If `dims` is a length-3 array, the array will have size `(dims...)`. 
"""
function amr_volume(snap; iv::Union{Int, Array, String} = 0, all = true,
                    dims::Union{Tuple, Int}=100, verbose=0)


    if typeof(dims) <: Tuple
        nx, ny, nz = dims
    else
        nx = ny = nz = dims
    end

    buffer = init_buffer(snap, iv, dims, 3)
    patches = snap["patches"]

    Size = snap["cartesian"]["size"]
    origin = snap["cartesian"]["origin"]

    x = range(origin[1], origin[1]+Size[1], length=nx)
    y = range(origin[2], origin[2]+Size[2], length=ny)
    z = range(origin[3], origin[3]+Size[3], length=nz)
    for iv in keys(buffer)
        for patch in patches

            x_extent = patch["x"]
            y_extent = patch["y"]
            z_extent = patch["z"]

            x_ids = findall(x_extent[1] .< x .< x_extent[end])
            y_ids = findall(y_extent[1] .< y .< y_extent[end])
            z_ids = findall(z_extent[1] .< z .< z_extent[end])

            itp = Itp.interpolate((x_extent, y_extent, z_extent), box(patch, iv=iv, verbose=verbose, all=true), Itp.Gridded(Itp.Linear()))
            @inbounds for iz in z_ids, iy in y_ids, ix in x_ids
                buffer[iv][ix, iy, iz] = itp(x[ix], y[iy], z[iz])
            end
        end
    end

    if length(keys(buffer)) == 1
        return collect(values(buffer))[1]
    else
        return buffer
    end
end

"""
    unigrid_plane(snap::Dict; iv::Union{Int, String}, x::Float, y::Float, z::Float)

Return a 2d array of joined patch data in a slice at `x/y/z`. All patches in given plane must have same shape, size,
and number of cells. Returns a `Dict` if `iv` is a collection, otherwise a `Matrix`.

#Arguments:
- `snap::Dict`, snapshot object

"""
function unigrid_plane(snap::Dict; x = nothing, y = nothing, z = nothing,
                      iv = 0, verbose=0, all=false)

    xyz = [x, y, z]
    ax = getindex((1, 2, 3), xyz .!= nothing)[1]

    patches = patches_in(snap, x=x, y=y, z=z)
    if length(patches) == 0
      throw(ErrorException(" no patches found in [$x, $y, $z]"))
    end

    verbose > 0 && println("number of patches: $(length(patches))")

    n1, n2 = 0, 0
    patchDict = Dict{Int, Tuple{NTuple{4,Int64}, Dict}}()
    for p in patches
        idxs = corner_indices(snap, p, dir=ax)
        if idxs[2] > n1
            n1 = idxs[2]
        end

        if idxs[4] > n2
            n2 = idxs[4]
        end

        patchDict[p["id"]] = (idxs, p)
    end
    datashp = [n1, n2]

    verbose >= 1 && println("data shape: $datashp") 


    buffer = init_buffer(snap, iv, datashp, 2)
    for iv in keys(buffer)
        patchkeys = keys(patchDict) |> collect
        Base.Threads.@threads for key in patchkeys
            idxs, patch = patchDict[key]
            im = plane(patch, x = x, y = y, z = z, iv = iv,
                             verbose=verbose, all=all)
            
            buffer[iv][idxs[1]:idxs[2], idxs[3]:idxs[4]] = im
        end
    end

    if length(keys(buffer)) == 1
        return collect(values(buffer))[1]
    else
        return buffer
    end

end



"""
    unigrid_volume(snap::Dict; iv::Union{Int, String}=0, all::Bool=false, verbose::Int=0)

Return a 3D array of joined patch data of quantity . All patches must have same shape, size, and
number of cells.

#Arguments:
- `snap::Dict`, snapshot

#Kwargs:
    - `iv::Union{String, Int, Array}`, what quantity(/ies) to buffer. if array, output will
          be a dictionairy with keys equal to the input variable names, and
          values as the 2d-buffers.
    - `all::Bool`, whether to include guard zones (if available)
"""
function unigrid_volume(snap; iv = 0,  all=false, verbose=0)

    dims = snap["cartesian"]["dims"]
    patches = snap["patches"]

    patchDict = Dict{Int, Tuple{NTuple{6, Int64}, Dict}}()

    nx, ny, nz = 0, 0, 0
    for p in patches
        idxs_xy = corner_indices(snap, p, dir=3)
        idxs_z = corner_indices(snap, p, dir=1)[3:4]

        idxs_xy[2] > nx ? nx = idxs_xy[2] : nothing
        idxs_xy[4] > ny ? ny = idxs_xy[4] : nothing
        idxs_z[2]  > nz ? nz = idxs_z[2]  : nothing

        idxs = tuple(idxs_xy..., idxs_z...)
        patchDict[p["id"]] = (idxs, p)
    end

    datashp = (nx, ny, nz)
    buffer = init_buffer(snap, iv, datashp, 3)
    verbose == 1 && println("volume shape $datashp")

    for iv in keys(buffer)
        for (idxs, patch) in values(patchDict)
            data = box(patch, iv=iv, all=all, verbose=verbose)

            buffer[iv][idxs[1]:idxs[2],
                       idxs[3]:idxs[4],
                       idxs[5]:idxs[6]] = data

        end
    end

    if length(keys(buffer)) == 1
        return collect(values(buffer))[1]
    else
        return buffer
    end
end

"""
    resample(d1::Array{Float, 1} where T, d2::Array{Float, 1},
             data::Array{Float, 2}, newdims::Union{Int, Tuple})

Resize 2d data array to shape defined by `newdims`. Data is resampled over axis values d1 and d2
using gridded bilinear interpolation.

Arguments:
--------------
- `d1::Array{Float, 1}`, axis values along dimension 1
- `d2::Array{Float, 1}`, axis values along dimension 2
- `data::Array{Float, 2}`, data to be resampled. Must have `shape == (length(d1), length(d2))`
- `newdims::Union{Int, Tuple}`, shape of resampled data. If `Int`, all dimensions will have
               same size. If `Tuple`, `size(newdims)` must be 2.

Returns:
--------------
- `Array{Float, 1}`, resampled axis 1 values
- `Array{Float, 1}`, resampled axis 2 values
- `Array{Float, 2}`, resampled data values
"""
function resample(d1, d2, data::Array{T, 2} where T, newdims)

    if typeof(newdims) <: Int
        n1, n2 = newdims, newdims
    else
        n1, n2 = newdims
    end

    nd1 = range(d1[1], d1[end], length=n1)
    nd2 = range(d2[1], d2[end], length=n2)

    new_data = similar(data, n1, n2)

    itp = Itp.interpolate((d1, d2), data, Itp.Gridded(Itp.Linear()))

    for i2 in eachindex(nd2)
        for i1 in eachindex(nd1)
            new_data[i1, i2] = itp(nd1[i1], nd2[i2])
        end
    end

    return nd1, nd2, new_data
end

"""
    resample(xs::Array{Float, 1} where T, ys::Array{Float, 1}, zs::Array{Float, 1},
            data::Array{Float, 3}, newdims::Union{Int, Tuple})

Resize 3d data array to shape defined by newdims. Data is resampled over axis values xs, ys, and zs
using gridded trilinear interpolation.

Arguments:
--------------
- `xs::Array{Float, 1}`, axis values along dimension 1
- `ys::Array{Float, 1}`, axis values along dimension 2
- `zs::Array{Float, 1}`, axis values along dimension 3
- `data::Array{Float, 3}`, data to be resampled. Must have `shape == (length(d1), length(d2), length(d3))`
- `newdims::Union{Int, Tuple}`, shape of resampled data. If `Int`, all dimensions will have
               same size. If `Tuple`, `size(newdims)` must be 3.

Returns:
--------------
- `Array{Float, 1}`, resampled axis 1 values
- `Array{Float, 1}`, resampled axis 2 values
- `Array{Float, 1}`, resampled axis 3 values
- `Array{Float, 3}`, resampled data values
"""
function resample(xs, ys, zs, data::Array{T, 3} where T, newdims)

    if typeof(newdims) <: Int
        nx = ny = nz = newdims
    else
        nx, ny, nz = newdims
    end

    nxs = range(xs[1], xs[end], length=nx)
    nys = range(ys[1], ys[end], length=ny)
    nzs = range(zs[1], zs[end], length=nz)

    new_data = similar(data, nx, ny, nz)

    itp = Itp.interpolate((xs, ys, zs), data, Itp.Gridded(Itp.Linear()))

    for iz in eachindex(nzs)
        for iy in eachindex(nys)
            for ix in eachindex(nxs)
                new_data[ix, iy, iz] = itp(nxs[ix], nys[iy], nzs[iz])
            end
        end
    end

    return nxs, nzs, nzs, new_data
end

