

function interpolate_3d(patch, data; x = 0.5, y = 0.5, z=0.5, verbose = 0)
    """
    Interpolate 3-dimensional data on to a point [x, y, z]

    Input:
    -----------------------------
    patch: dict, patch from a given snapshot
    data: array of floats, 3d-data for a given variable
    x, y, z: float, point at which to interpolate

    Output:
    - value of data at [x, y, z]
    """

    nx, ny, nz = size(data)
    xs, ys, zs = patch["x"], patch["y"], patch["z"]

    # find the points below the integration point
    ix, iy, iz = 0, 0, 0
    try
        ix = findall(xs .< x)[end]
        iy = findall(ys .< y)[end]
        iz = findall(zs .< z)[end]
    catch
        nothing
    end

    # weights
    xd = 1.0/(xs[ix+1] - xs[ix])*(x - xs[ix])
    yd = 1.0/(ys[iy+1] - ys[iy])*(y - ys[iy])
    zd = 1.0/(zs[iz+1] - zs[iz])*(z - zs[iz])

    c000 = data[ix  , iy  , iz]
    c100 = data[ix+1, iy  , iz]
    c110 = data[ix+1, iy+1, iz]
    c111 = data[ix+1, iy+1, iz+1]
    c101 = data[ix+1, iy  , iz+1]
    c010 = data[ix  , iy+1, iz]
    c001 = data[ix  , iy  , iz+1]
    c011 = data[ix  , iy+1, iz+1]

    c00 = c000*(1 - xd) + c100*xd
    c01 = c001*(1 - xd) + c101*xd
    c10 = c010*(1 - xd) + c110*xd
    c11 = c011*(1 - xd) + c111*xd

    c0 = c00*(1 - yd) + c10*yd
    c1 = c01*(1 - yd) + c11*yd

    itpPoint = c0*(1 - zd) + c1*zd
    return itpPoint

end

function interpolate_2d(patch::Dict, data::AbstractArray; ax = 1, p1 = 0.5, p2 = 0.5, verbose = 0)

    n1, n2 = size(data)
    d1, d2 = getindex(("x", "y", "z"), [1, 2, 3] .!= ax)
    d1, d2 = patch[d1], patch[d2]

    # find the points below the integration point
    i1, i2 = 0, 0
    try
        i1 = findall(d1 .< p1)[end]
        i2 = findall(d2 .< p2)[end]
    catch
        nothing
    end

    # coefficients
    c11 = data[i1, i2]
    c21 = data[i1+1, i2]
    c12 = data[i1, i2+1]
    c22 = data[i1+1, i2+1]

    w1 = 1.0/(d1[i1+1] - d1[i1])*(p1 - d1[i1])
    w2 = 1.0/(d2[i2+1] - d2[i2])*(p2 - d2[i2])

    f1 = (d1[i1+1] - p1)/(d1[i1+1] - d1[i1])*c11 + w1*c21
    f2 = (d1[i1+1] - p1)/(d1[i1+1] - d1[i1])*c12 + w1*c22

    f = (d2[i2+1] - p2)/(d2[i2+1] - d2[i2])*f1 .+
        (p2 - d2[i2])/(d2[i2+1] - d2[i2])*f2

    return f
end

"""
    interpolate(patch::Dict; iv::Union{Int, String}, x::Float, y::Float, z::Float,
                Log::Bool, all::Bool)

Find the plane values of a patch at a given slice x/y/z by interpolating
neighbouring planes.

# Arguments:
- `patch::Dict`, patch object from a snapshot

# Kwargs:
- `iv::Union{Int, String}`: quantity to get data of, default `0`
- `(x, y, z)::Float`: position at which to slice, default `nothing`
- `Log::Bool`: whether to log the data, default `false`
- `all::Bool`: whether to include guard zones, default `false`

# Returns:
- 2d array of float32, interpolated values at position x/y/z

"""
function interpolate(patch; iv = 0, x = nothing, y = nothing, z = nothing,
                     Log=false, all=false)

    if patch["guard_zones"]
        li = patch["li"]  # lower inner
        ui = patch["ui"]  # upper inner
    else
        li = ones(Int, 3) # using Static??
        ui = patch["n"]
    end

    xyz = [x, y, z]
    dir = getindex((1, 2, 3), xyz .!= nothing)[1]
    dirStr = getindex(("x", "y", "z"), xyz .!= nothing)[1]
    ax = xyz[xyz .!= nothing][1]

    n = patch["n"][dir]

    i = findall(patch[dirStr] .< ax)[end]

    w = 1.0/(patch[dirStr][i+1] - patch[dirStr][i]) *
        (ax - patch[dirStr][i])

    data = patch["var"](iv)

    if iv == 0 || iv == 4 || iv == "d" || iv == "e"
        Log = true
    end

    if Log
        data = log.(data)
    end

    if !all
        if dir == 1
            data = (1 - w) .* data[i,   li[2]:ui[2], li[3]:ui[3]] .+
                        w  .* data[i+1, li[2]:ui[2], li[3]:ui[3]]
        elseif dir == 2
            data = (1 - w) .* data[li[1]:ui[1],i,li[3]:ui[3]] .+
                        w  .* data[li[1]:ui[1],i+1,li[3]:ui[3]]
            data = data'
        else
            data = (1 - w) .* data[li[1]:ui[1],li[2]:ui[2],i] .+
                        w  .* data[li[1]:ui[1],li[2]:ui[2],i+1]
        end
    else
        if dir == 1
            data = (1 - w) .* data[i,   :, :] .+
                        w  .* data[i+1, :, :]
        elseif dir == 2
            data = (1 - w) .* data[:,i,:] .+
                        w  .* data[:,i+1,:]
            data = data'
        else
            data = (1 - w) .* data[:,:,i] .+
                        w  .* data[:,:,i+1]
        end
    end

    Log ? exp.(data) : data
end
