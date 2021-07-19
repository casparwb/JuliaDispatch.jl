

using JuliaDispatch.Select, Unitful
using Interpolations
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
        datashp = reverse(dims) # or array of points in each dimension
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
             Log::Bool, dims::Union{Int, Tuple}, verbose=0)

Return a 2D array containing interpolated data of quantity `iv` in a slice `x/y/z` from all patches in a given snapshot. If `dims` is an
`Int` the resulting array will be of size `(dims, dims)`. If `dims` is a length-2 array, the array will have size `(dims...)`. 
"""
function amr_plane(snap; iv = 0, x = nothing, y = nothing, z = nothing,
                   Log = false, dims::Union{Int, Tuple}=100, all=false, span=nothing, verbose=0)


    xyz = [x, y, z]
    axIdx = getindex((1, 2, 3), xyz .!= nothing)[1]
    ax = xyz[axIdx]
    dir1, dir2 = getindex((1, 2, 3), xyz .== nothing)
    dir1s, dir2s = getindex(("x", "y", "z"), xyz .== nothing)


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

    axis1 = range(origin[1], origin[1]+Size[1], length=n1) # axis in plane axis 1
    axis2 = range(origin[2], origin[2]+Size[2], length=n2) # axis in plane axis 2
    if !isnothing(span)
        axis1_span, axis2_span = span

        snapshot_span = ((axis1[1], axis1[end]), (axis2[1], axis2[end]))

        # check if given span extends beyond boundary and if corresponding axis is periodic
        axis1_over, axis1_under = false, false
        if (axis1_span[1] < snapshot_span[1][1])  
            axis1_under = true
            verbose == 2 && println("given span extends under axis 1 boundary")
        elseif (axis1_span[2] > snapshot_span[1][2])
            axis1_over = true
            verbose == 2 && println("given span extends over axis 1 boundary")

        end

        any((axis1_over, axis1_under)) && iszero(snap["periodic"][dir1]) && throw(ArgumentError("axis $(dir1) is not periodic"))

        axis2_over, axis2_under = false, false
        if (axis2_span[1] < snapshot_span[2][1])  
            axis2_under = true
            verbose == 2 && println("given span extends under axis 2 boundary")
        elseif (axis2_span[2] > snapshot_span[2][2])
            axis2_over = true
            verbose == 2 && println("given span extends over axis 2 boundary")
        end

        any((axis2_over, axis2_under)) && iszero(snap["periodic"][dir2]) && throw(ArgumentError("axis $(dir2) is not periodic"))

        sub_axis1 = range(axis1_span[1], axis1_span[2], length=n1) |> collect # axis in plane axis 1
        sub_axis2 = range(axis2_span[1], axis2_span[2], length=n2) |> collect # axis in plane axis 2

        if any((axis1_under, axis1_over, axis2_under, axis2_over)) 
            
            if axis1_under
                ids = findall(sub_axis1 .< snapshot_span[1][1])
                stopidx =  findall(axis1 .>= axis1[1] - sub_axis1[1])[1]
                sub_axis1[ids] = range(axis1[end - stopidx], stop=axis1[end], length=length(ids))
            elseif axis1_over
                ids = findall(sub_axis1 .> snapshot_span[1][2])
                stopidx = findall(axis1 .<= sub_axis1[end] - axis1[end])[end]
                sub_axis1[ids] = range(axis1[1], stop=axis1[stopidx], length=length(ids))
            end

            if axis2_under
                ids = findall(sub_axis2 .< snapshot_span[2][1])
                stopidx =  findall(axis2 .>= axis2[1] - sub_axis2[1])[1]
                sub_axis2[ids] = range(axis2[end - stopidx], stop=axis2[end], length=length(ids))
            elseif axis2_over
                ids = findall(sub_axis2 .> snapshot_span[2][2])
                stopidx = findall(axis2 .<= sub_axis2[end] - axis2[end])[end]
                sub_axis2[ids] = range(axis2[1], stop=axis2[stopidx], length=length(ids))
            end

        end

        axis1 = sub_axis1
        axis2 = sub_axis2
    end

    interpx(itp, y, z) = itp(ax, y, z)
    interpy(itp, x, z) = itp(x, ax, z)
    interpz(itp, x, y) = itp(x, y, ax)
    interp = [interpx, interpy, interpz][axIdx] 


    for iv in keys(buffer)
        Base.Threads.@threads for patchID in eachindex(patches)
            patch = patches[patchID]

            if patch["guard_zones"] && !all
                li = patch["li"] .- 1  # lower inner
                ui = patch["ui"] .+ 1  # upper inner
            elseif patch["guard_zones"] && all
                li = ones(Int, 3)
                ui = patch["gn"]
            else
                li = ones(Int, 3)
                ui = patch["n"]
            end

            axis1extent = patch[dir1s][li[dir1]:ui[dir1]]
            axis2extent = patch[dir2s][li[dir2]:ui[dir2]]

            axis1_indices = findall(axis1extent[1] .< axis1 .< axis1extent[end]) # indices extended by patch in plane axis 1
            axis2_indices = findall(axis2extent[1] .< axis2 .< axis2extent[end]) # indices extended by patch in plane axis 2

            itp = LinearInterpolation((patch["x"], patch["y"], patch["z"]), box(patch, iv=iv, all=true), extrapolation_bc=Throw())
            for i2 ∈ axis2_indices, i1 ∈ axis1_indices
                buffer[iv][i2, i1] = interp(itp, axis1[i1], axis2[i2]) # interpolate and insert into buffer
            end
        end
    end

    if length(keys(buffer)) == 1
        return axis1, axis2, collect(values(buffer))[1]
    else
        return axis1, axis2, buffer
    end

end


"""
    amr_volume(snap::Dict; iv::Union{Int, AbstractArray, String}, all::Bool,
            dims::Union{Tuple, Int}, verbose::int)

Return a 3D array containing interpolated data of quantity `iv` from all patches in a given snapshot. If `dims` is an
`Int` the resulting array will be of size `(dims, dims, dims)`. If `dims` is a length-3 array, the array will have size `(dims...)`. 
"""
function amr_volume(snap; iv::Union{Int, Array, String} = 0, all = false,
                    dims::Union{Tuple, Int}=100, verbose=0, span=nothing)


    if typeof(dims) <: Tuple
        nx, ny, nz = dims
    else
        nx = ny = nz = dims
    end

    buffer = init_buffer(snap, iv, dims, 3)
    
    Size = snap["cartesian"]["size"]
    origin = snap["cartesian"]["origin"]
    
    if isnothing(span)
        x = range(origin[1], origin[1]+Size[1], length=nx)
        y = range(origin[2], origin[2]+Size[2], length=ny)
        z = range(origin[3], origin[3]+Size[3], length=nz)
        patches = patches_in(snap, span)
        # patches = snap["patches"]
    else
        x_span, y_span, z_span = span
        x = range(x_span[1], x_span[2], length=nx)
        y = range(y_span[1], y_span[2], length=ny)
        z = range(z_span[1], z_span[2], length=nz)
        patches = snap["patches"]
    end

    for iv in keys(buffer)
        Base.Threads.@threads for patch in patches

            if patch["guard_zones"] && !all
                li = patch["li"] .- 1  # lower inner
                ui = patch["ui"] .+ 1  # upper inner
            elseif patch["guard_zones"] && all
                li = ones(Int, 3)
                ui = patch["gn"]
            else
                li = ones(Int, 3)
                ui = patch["n"]
            end

            x_extent = patch["x"][li[1]:ui[1]]
            y_extent = patch["y"][li[2]:ui[2]]
            z_extent = patch["z"][li[3]:ui[3]]

            x_ids = findall(x_extent[1] .< x .< x_extent[end])
            y_ids = findall(y_extent[1] .< y .< y_extent[end])
            z_ids = findall(z_extent[1] .< z .< z_extent[end])

            itp = LinearInterpolation((patch["x"], patch["y"], patch["z"]), box(patch, iv=iv, verbose=verbose, all=true), extrapolation_bc=Throw())
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

"""
function unigrid_plane(snap::Dict; x = nothing, y = nothing, z = nothing,
                      iv = 0, verbose=0, all=false, span = nothing)

    xyz = [x, y, z]
    ids = [1, 2, 3]
    ax = getindex(ids, xyz .!= nothing)[1]
    deleteat!(ids, ax)

    patches = patches_in(snap, x=x, y=y, z=z)
    if length(patches) == 0
      throw(ErrorException(" no patches found in [$x, $y, $z]"))
    end

    verbose > 0 && println("number of patches: $(length(patches))")

    n1, n2 = 0, 0
    patchDict = Dict{Int, Tuple{NTuple{4,Int64}, Dict}}()
    for p in patches
        idxs = corner_indices(snap, p, dir=ax)

        n1 = max(n1, idxs[2])
        n2 = max(n2, idxs[4])

        patchDict[p["id"]] = (idxs, p)
    end
    datashp = [n1, n2]

    verbose >= 1 && println("data shape: $datashp") 

    if !isnothing(span)
        span = (span[1] .* 1.0, span[2] .* 1.0)
        snapshot_span = [(snap["cartesian"]["origin"][i], snap["cartesian"]["origin"][i] + snap["cartesian"]["size"][i]) for i = 1:3]

        # check if given span extends beyond boundary and if corresponding axis is periodic
        axis1_span, axis2_span = span
        if ((snapshot_span[ids[1]][1] > axis1_span[1]) || (snapshot_span[ids[1]][2] < axis1_span[2]))    
            iszero(snap["periodic"][ids[1]]) && throw(ArgumentError("axis $(ids[1]) is not periodic"))
        end

        if ((snapshot_span[ids[2]][1] > axis2_span[1]) || (snapshot_span[ids[2]][2] < axis2_span[2]))
            iszero(snap["periodic"][ids[2]]) && throw(ArgumentError("axis $(ids[2]) is not periodic"))
        end

        span = Array([span...])
        span = tuple(insert!(span, ax, snapshot_span[ax])...)

        max_ids = [1, 1]
        min_ids = [n1, n1]

        for patch in patches_in(snap, span)
            idxs = corner_indices(snap, patch, dir=ax)

            max_ids[1] = max(max_ids[1], idxs[2])
            max_ids[2] = max(max_ids[2], idxs[4])

            min_ids[1] = min(min_ids[1], idxs[1])
            min_ids[2] = min(min_ids[2], idxs[3])
        end

    else
        max_ids = [n1, n2]
        min_ids = [1, 1]

    end


    buffer = init_buffer(snap, iv, datashp, 2)
    for iv in keys(buffer)
        patchkeys = keys(patchDict) |> collect
        Base.Threads.@threads for key in patchkeys
            idxs, patch = patchDict[key]
            im = plane(patch, x = x, y = y, z = z, iv = iv,
                             verbose=verbose, all=all)
            
            buffer[iv][idxs[3]:idxs[4], idxs[1]:idxs[2]] = im'
        end
    end

    if length(keys(buffer)) == 1
        return collect(values(buffer))[1][min_ids[2]:max_ids[2], min_ids[1]:max_ids[1]]
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

