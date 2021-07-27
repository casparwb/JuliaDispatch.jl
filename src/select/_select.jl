using Interpolations

using StaticArrays
"""
    values_along(snap::Dict, point::Array=[0.5, 0.5, 0.5]; dir::Int=1, iv::Union{String, Int}=0, var=nothing,
                    verbose::Int = 0, all::Bool = false)


Return `s, f(s)` with `s` the coordinates and `f` the values in the `iv`
slot of data, taken along the direction v -- so far restricted
to axis values
"""
function values_along(snap, point=[0.5, 0.5, 0.5];
                    dir=1, iv=0, var=nothing, verbose = 0, all = false)
                    

    patches = patches_along(snap, point, dir = dir, verbose = verbose)

    if isnothing(patches)
        throw(ErrorException("no patches found in [$(point[1]), $(point[2]), $(point[3])]"))
    end

    ss = Array{Float32}([])
    ff = Array{Float32}([])

    @inbounds for p in patches
        ss1, ff1 = values_in(p, point, dir=dir, iv=iv, var=var, all=all, verbose=verbose)
        append!(ss, ss1)
        append!(ff, ff1)
    end

    return ss, ff

end

""" 
    patches_in(snap::Dict; x::Number=nothing, y::Number=nothing, z::Number=nothing)

Return a Vector of patches that lie within a given plane. 
"""
function patches_in(snap; x = nothing, y = nothing, z = nothing)

    if !isnothing(x)
        patches = [p for p in snap["patches"]
                if (x >= p["extent"][3, 1] && x <= p["extent"][3, 2])]
    end

    if !isnothing(y)
        patches = [p for p in snap["patches"]
                if (y >= p["extent"][1, 1] && y <= p["extent"][1, 2])]
    end

    if !isnothing(z)
        patches = [p for p in snap["patches"]
                if (z >= p["extent"][2, 1] && z <= p["extent"][2, 2])]
    end

    return patches
end

""" 
    patches_in(snap::Dict, span::Tuple)

Return a Vector of patches that lie within a domain spanned by `span`.
"""
function patches_in(snap, span)

    xspan, yspan, zspan = span
    patches = Dict[]
    # println(xspan, " ", yspan, " ", zspan)
    for p in snap["patches"]
        if ((xspan[1] <= p["extent"][3, 1] && xspan[2] >= p["extent"][3, 2]) &&
            (yspan[1] <= p["extent"][1, 1] && yspan[2] >= p["extent"][1, 2]) &&
            (zspan[1] <= p["extent"][2, 1] && zspan[2] >= p["extent"][2, 2]))
            push!(patches, p)
        end
    end

    return patches
end


"""
    values_in(snap::Dict, point::Array=[0.5, 0.5, 0.5]; dir::Int=1, iv::Union{String, Int}=0, var=nothing,
    verbose::Int = 0, all::Bool = false)


Return `s, f(s)` with `s` the coordinates and `f` the values in the `iv`
slot of data, taken along the direction v -- so far restricted
to axis values
"""
function values_in(p, point = [0.5,0.5,0.5];
                    dir = 1, iv = 0, i4 = 1, var = nothing, verbose = 0, all = false)


    ss, ff = [], []

    ii, w = indices_and_weights(p, point, iv)
    data = p["var"](iv, i4=i4, all=all, verbose=verbose)

    # +1 because of julia 1-indexing
    ione = (0, 1)[(p["gn"][1] > 1) + 1]
    jone = (0, 1)[(p["gn"][2] > 1) + 1]
    kone = (0, 1)[(p["gn"][3] > 1) + 1]


    if !p["guard_zones"] || !all
        m = @SVector ones(Int32, 3)
        n = p["n"]
    elseif all
        m = @SVector ones(Int32, 3)
        n = p["gn"]
    end


    if typeof(var) <: String
        iv = p["idx"]["dict"][var]
    end

    if dir == 1
        j = ii[2]
        k = ii[3]

    for i ∈ m[1]:m[1]+n[1]-1
        f1 = data[i, j, k     ]*(1 - w[2]) + data[i, j+jone, k     ]*w[2]
        f2 = data[i, j, k+kone]*(1 - w[2]) + data[i, j+jone, k+kone]*w[2]
        f  = f1*(1 - w[3]) + f2*w[3]

        if "p1" in keys(p)
            iv == p["idx"]["p1"] ? push!(ss, p["xs"][i]) : push!(ss, p["x"][i])
        else
            push!(ss, p["x"][i])
        end
        push!(ff, f)
    end

    elseif dir == 2
    i = ii[1]
    k = ii[3]

    for j = m[2]:m[2]+n[2]-1
        f1 = data[i, j, k     ]*(1-w[1]) + data[i+ione, j, k     ]*w[1]
        f2 = data[i, j, k+kone]*(1-w[1]) + data[i+ione, j, k+kone]*w[1]
        f  = f1*(1-w[3]) + f2*w[3]
        if "p2" in keys(p)
            iv == p["idx"]["p2"] ? push!(ss, p["ys"][j]) : push!(ss, p["y"][j])
        else
            push!(ss, p["y"][j])
        end
        push!(ff, f)
    end

    else
        i = ii[1]
        j = ii[2]

        for k = m[3]:m[3]+n[3]-1
            f1 = data[i, j     , k]*(1-w[1]) + data[i+ione, j     , k]*w[1]
            f2 = data[i, j+jone, k]*(1-w[1]) + data[i+ione, j+jone, k]*w[1]
            f  = f1*(1-w[2]) + f2*w[2]
        if "p3" in keys(p)
            iv == p["idx"]["p3"] ? push!(ss, p["zs"][k]) : push!(ss, p["z"][k])
        else
            push!(ss, p["z"][k])
        end
        push!(ff, f)
    end

    end

    return ss, ff

end

"""
    indices_and_weights(p::Dict, point = [0.5, 0.5, 0.5], iv = 0)

Return indices and interpolation weights for a point in a patch p.
"""
function indices_and_weights(p::Dict, point = [0.5, 0.5, 0.5], iv = 0)

    x0 = p["x"][1]
    y0 = p["y"][1]
    z0 = p["z"][1]

    if "p1" in keys(p["idx"])
        iv == p["idx"]["p1"] ? x0 = p["xs"][1] : nothing
        iv == p["idx"]["p2"] ? y0 = p["ys"][1] : nothing
        iv == p["idx"]["p3"] ? z0 = p["zs"][1] : nothing
    end
    if "b1" in keys(p["idx"])
        iv == p["idx"]["b1"] ? x0 = p["xs"][1] : nothing
        iv == p["idx"]["b2"] ? y0 = p["ys"][1] : nothing
        iv == p["idx"]["b3"] ? z0 = p["zs"][1] : nothing
    end

    corner = SVector{3, Float32}(x0, y0, z0)
    weights = SVector{3, Float32}(point .- corner) ./ p["ds"]
    indices = SVector{3, Int32}(round.(weights))

    # indices = [max(0, min(p["n"][i] - 2, indices[2])) for i = 1:3]
    indices = SVector{3, Int32}([max(0, min(p["n"][i] - 2, indices[2])) for i = 1:3] .+ 1)

    # indices .+= 1

    return indices, weights

end


"""
    patches_along(snap, point=[0.5, 0.5, 0.5]; dir = 1, verbose=0)

Get the patches along a given direction through a point
"""
function patches_along(snap, point=[0.5, 0.5, 0.5]; dir = 1, verbose=0)

    if haskey(snap, "patches")#"patches" in keys(snap) 
        patches = snap["patches"]
    end

    pt = copy(point)

    p = patch_at(patches, pt)

    if !isnothing(p)
        p1 = p
        while !isnothing(p)
            pt[dir] = p["position"][dir] - p["size"][dir]*0.6
            p = patch_at(patches, pt)
            if !isnothing(p)
                p1 = p
            end
        end

        p = p1
        out = []
        while !isnothing(p)
            push!(out, p)
            pt[dir] = p["position"][dir] + p["size"][dir]*0.6
            p = patch_at(patches, pt)
        end
    else
        out = nothing
    end

    return out

end

"""
    patch_at(pp, point=[0.5, 0.5, 0.5]; verbose=0)

Find the patch that contains a given point.
"""
function patch_at(patches, point=[0.5, 0.5, 0.5]; verbose=0)

    level = -1
    p1 = nothing
    @inbounds for p in patches
        if is_inside(p, point, verbose=verbose)
            p["level"] > level ? p1 = p : level = p["level"]
        end
    end

    return p1

end

"""
    is_inside(p::Dict, point::Array)

Return true/false depending on if `point` is inside `patch`.
"""
function is_inside(p, point; verbose=0)

    verbose == 1 && println("""checking in $point is inside patch $(patch["id"])""")
    left = point .>= p["llc_cart"]
    right = point .<= p["llc_cart"] + p["size"]

    all(left) && all(right)

end

"""
    corner_indices(snap, patch; dir=1)

Get the corner indices of a given patch in a given direction.
"""
function corner_indices(snap, patch; dir=1)


    i = (patch["position"] .- patch["size"]/2 .- snap["cartesian"]["origin"])./patch["ds"]
    i = [Int(round(k + 0.5)) for k in i]
    n = patch["n"][:]

    if dir >= 1 && dir < 4
        deleteat!(i, dir)
        deleteat!(n, dir)
        return i[1]+1, i[1]+n[1], i[2]+1, i[2]+n[2]
    else
        return i[1], i[1]+n[1], i[2], i[2]+n[2], i[3], i[3]+n[3]
    end
end

"""
    corner_indices(snap, patch)

Get the corner indices of a given patch.
"""
function corner_indices(snap, patch)


    i = (patch["position"] .- patch["size"]/2 .- snap["cartesian"]["origin"])./patch["ds"]
    i = [Int(round(k + 0.5)) for k in i]
    n = patch["n"][:]

    # if dir >= 1 && dir < 4
    #     deleteat!(i, dir)
    #     deleteat!(n, dir)
    #     return i[1]+1, i[1]+n[1], i[2]+1, i[2]+n[2]
    # else
    # end
    return i[1]+1, i[1]+n[1], i[2]+1, i[2]+n[2], i[3]+1, i[3]+n[3]
end


"""
box(patch::Dict; iv::Union{String, Int}, all::Bool, verbose::Int)

Return 3D-array of data of quantity `iv` from patch. Shorthand for `patch["var"](iv, verbose=verbose, all=all)`.

# Examples
```
julia> density_box = box(patch, iv="d")
```
"""
function box(patch; iv = 0, all=false, verbose = 0)
    f = patch["var"](iv, verbose=verbose, all=all)
end


"""
    plane(patch::Dict; iv::Union{Int, String}, x::Float, y::Float, z::Float,
                Log::Bool, all::Bool)

Find the plane values of a patch at a given slice x/y/z by interpolating
neighbouring planes.

# Examples
```
julia> density_plane = plane(patch, iv="d", z=0.0)
```
"""
function plane(patch; iv = 0, x = nothing, y = nothing, z = nothing,
                     Log=false, all=false, verbose=0)

    if patch["guard_zones"] && !all
        li = patch["li"]  # lower inner
        ui = patch["ui"]  # upper inner
    elseif patch["guard_zones"] && all
        li = ones(Int, 3)
        ui = patch["gn"]
    else
        li = ones(Int, 3)
        ui = patch["n"]
    end

    data = box(patch, iv=iv, all=true)
    if any(map(isone, size(data)))
        idx = findall(isone, size(data))[1]
        return selectdim(data, idx, 1)
    end

    itp = LinearInterpolation((patch["x"], patch["y"], patch["z"]), 
                               data, 
                               extrapolation_bc = Throw())

    if !isnothing(x)
        verbose == 1 && println("plane at x = $x")
        return itp(x, patch["y"], patch["z"])[li[2]:ui[2], li[3]:ui[3]] # interpolate in x-slice

    elseif !isnothing(y)
        verbose == 1 && println("plane at y = $y")
        return itp(patch["x"], y, patch["z"])[li[1]:ui[1], li[3]:ui[3]] # interpolate in y-slice

    elseif !isnothing(z)
        verbose == 1 && println("plane at z = $z")
        return itp(patch["x"], patch["y"], z)[li[1]:ui[1], li[2]:ui[2]] # interpolate in z-slice
    end
end


