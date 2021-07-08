

using JuliaDispatch.Select, JuliaDispatch.Interpolations, Unitful
import Interpolations
const Itp = Interpolations

"""
init_buffer(snap::Dict, iv::Union{Int, String}, dims::Union{Tuple, AbstractArray},
            ndims::Int)

Function for initializing a buffer array for storing the data.

Arguments:
---------------
    - snap: Dict, snapshot dictionairy

Kwargs:
---------------------------
    - iv:    String/Int/Collection of String/Ints, what quantity(/ies) to extract, default 0
    - dims:  Int/Tuple/Array, data size in each dimension. If Int,
             all dimensions will have same length. If Tuple/Array,
             must have length(dims) == ndims.
    - ndims: number of dimensions

Output:
---------------------------
    - buffer: If iv is collection, output will be dictionairy with keys
              equal to the input iv names, and values as undefined
              n-dimensional arrays of type Float32.
              If iv is int/string, output is undefined ndims-dimensional array of Float32.
"""
function init_buffer(snap, iv, dims, ndims)

    datashp = nothing
    if typeof(dims) <: Int
        datashp = repeat([dims], ndims)
    else
        datashp = dims
    end

    buffer = Dict{Union{String, Int}, Array{Quantity{Float32}, ndims}}()
    if typeof(iv) <: String && iv == "all"
        ivs = [k for (k, iv) in snap["idx"]
               if (typeof(iv) <: Int && iv > 0)]

        for iv in ivs
            buffer[iv] = Array{Float32, ndims}(undef, datashp...)
        end
    elseif typeof(iv) <: Array
        ivs = iv
        for iv_ in iv
            buffer[iv_] = Array{Float32, ndims}(undef, datashp...)
        end
    else
        ivs = [iv]
        buffer[iv] = Array{Float32, ndims}(undef, datashp...)
    end

    return buffer

end

"""
    amr_plane(snap::Dict; iv::Union{Int, String}, x::Float, y::Float, z::Float,
             Log::Bool, dims::Union{Int, Tuple})

Return an interpolated 2d buffer of mesh-refined patch data with at a slice
x/y/z, with dimensions dims.

# Arguments:

- `snap::Dict`, snapshot object

# Kwargs:

- iv:      String/Int/Collection of String/Ints, what quantity(/ies) to extract, default 0
- x, y, z: Float, position at which to slice, default nothing
- all:     Bool, whether to include guard zones, default false
- Log:     Bool, whether to log the data, default false
- dims:    Int/Tuple, data size in each dimension. If Int,
            both1 dimensions will have same length. If Tuple/Array,
            must have length(dims) 2. Default 100

# Returns:

- Dictionairy of Array{Float32, 2} if iv is a collection,
    or Array{Float32, 2} if iv is Int/String
"""
function amr_plane(snap; iv = 0, x = nothing, y = nothing, z = nothing,
                   Log = false, dims::Union{Int, Tuple}=100)


    xyz = [x, y, z]
    ax = getindex((1, 2, 3), xyz .!= nothing)[1]
    dir1, dir2 = getindex((1, 2, 3), xyz .== nothing)
    dir1s, dir2s = getindex(("x", "y", "z"), xyz .== nothing)

    e1s, e2s = nothing, nothing
    if ax == 1 || ax == 2
        e1s = [1, 2]
        e2s = [3, 4]
    elseif ax == 3
        e1s = [3, 4]
        e2s = [1, 2]
    end

    Size = copy(snap["cartesian"]["size"])
    origin = copy(snap["cartesian"]["origin"])

    deleteat!(Size, ax)
    deleteat!(origin, ax)

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

    d1 = range(origin[1], origin[1]+Size[1], length=n1)
    d2 = range(origin[2], origin[2]+Size[2], length=n2)

    for iv in keys(buffer)
        for patch in patches
            # interpolate between planes
            data = interpolate(patch, iv=iv, x=x, y=y, z=z, all=true)
            e = patch["extent"]
            e1 = e[ax, e1s]
            e2 = e[ax, e2s]

            d1s = findall(e1[1] .<= d1 .<= e1[2])
            d2s = findall(e2[1] .<= d2 .<= e2[2])

            itp = Itp.interpolate((patch[dir1s], patch[dir2s]), data, Itp.Gridded(Itp.Linear()))
            @inbounds for i1 in d1s, i2 in d2s
                buffer[iv][i1, i2] = itp(d1[i1], d2[i2])
            end
        end
        buffer[iv] = buffer[iv]'
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

Return an interpolated 3d array from all mesh-refined patch data, with size
defined by dims.

Arguments:
-----------
    - snap: Dict, snapshot object

Kwargs:
-----------
    - iv:      String/Int/Collection of String/Ints, what quantity(ies) to extract, default 0
    - dims:    Int/Tuple, data size in each dimension. If Int,
               all dimensions will have same length. If Tuple/Array,
               must have length(dims) == 3. Default 100.

Returns:
-----------
    - Dictionairy of Array{Float32, 3} if iv is a collection,
      or Array{Float32, 3} if iv is Int/String
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
            data = box(patch, iv=iv, verbose=verbose, all=true)

            e = patch["extent"]
            ex = e[3,1:2]
            ey = e[3,3:4]
            ez = e[2,1:2]

            xs = findall(ex[1] .<= x .<= ex[2])
            ys = findall(ey[1] .<= y .<= ey[2])
            zs = findall(ez[1] .<= z .<= ez[2])

            itp = Itp.interpolate((patch["x"], patch["y"], patch["z"]), data, Itp.Gridded(Itp.Linear()))
            @inbounds for ix in xs, iy in ys, iz in zs
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

Return a 2d array of joined patch data. All patches in given plane must have same shape, size,
and number of cells.

Arguments:
-----------
    - snap: Dict, snapshot object

Kwargs:
-----------
    - iv:      String/Int/Collection of String/Ints, what quantity(ies) to extract, default 0
    - x, y, z: Float, at what position to slice, default nothing
    - all:     Bool, whether to include guard zones (if available)

Returns:
-----------
    - Dictionairy of Array{Float32, 2} if iv is a collection,
      or Array{Float32, 2} if iv is Int/String
"""
function unigrid_plane(snap::Dict; x = nothing, y = nothing, z = nothing,
                      iv = 0, verbose=0, all=false)

    xyz = [x, y, z]
    ax = getindex((1, 2, 3), xyz .!= nothing)[1]
    println(ax)
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
        #buffer[iv] = (buffer[iv])u"m/s"
        for (idxs, patch) in values(patchDict)
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
unigrid_volume(snap::Dict; iv::Union{Int, String}, all:Bool, verbose::Int)

Return a 3D array of joined patch data. All patches must have same shape, size, and
number of cells.

Arguments:
-----------
    - snap: dict, snapshot

Kwargs:
-----------
    - iv: string/int/array, what quantity(/ies) to buffer. if array, output will
          be a dictionairy with keys equal to the input variable names, and
          values as the 2d-buffers.
    - all: bool, whether to include guard zones (if available)

Output:
-----------
    - Dictionairy of Array{Float32, 3} if iv is a collection,
      or Array{Float32, 3} if iv is Int/String
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
    verbose == 1 ? println("volume shape $datashp") : nothing

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

    itp = Itp.interpolate((x, y, z), data, Itp.Gridded(Itp.Linear()))

    for iz in eachindex(nzs)
        for iy in eachindex(nys)
            for ix in eachindex(nxs)
                new_data[ix, iy, iz] = itp(nxs[ix], nys[iy], nzs[iz])
            end
        end
    end

    return nxs, nzs, nzs, new_data
end

