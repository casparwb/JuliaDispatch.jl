module Select

    export values_along, values_in, indices_and_weights, patches_along, 
           patch_at, is_inside, corner_indices, patches_in, box, plane

    using StaticArrays
    function values_along(pp, point=[0.5, 0.5, 0.5];
                        dir=1, iv=0, var=nothing, verbose = 0, all = 0)
        """
        Return s,f(s) with s the coordinates and f the values in the iv
        slot of data, taken along the direction v -- so far restricted
        to axis values
        """

        "patches" in keys(pp) ? pp = pp["patches"] : nothing

        patches = patches_along(pp, point, dir = dir, verbose = verbose)

        if patches == nothing
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

    """ Returns which patches are in a given plane """
    function patches_in(snap; x = nothing, y = nothing, z = nothing)
        patches = snap["patches"]

        if x != nothing
            patches = [p for p in patches
                    if (x >= p["extent"][3, 1] && x < p["extent"][3, 2])]
        end

        if y != nothing
            patches = [p for p in patches
                    if (y >= p["extent"][1, 1] && y < p["extent"][1, 2])]
        end

        if z != nothing
            patches = [p for p in patches
                    if (z >= p["extent"][2, 1] && z < p["extent"][2, 2])]
        end

        return patches
    end

    function values_in(p, point = [0.5,0.5,0.5];
                        dir = 1, iv = 0, i4 = 1, var = nothing, verbose = 0, all = 0)
        """
        Return s, f(s) with s the coordinates and f the values in the iv
        slot of data, taken along the direction v -- so far restricted to
        axis values
        """


        ss, ff = [], []

        ii, w = indices_and_weights(p, point, iv)
        data = p["var"](iv, i4=i4, all=all, verbose=verbose)

        # ijk = p["gn"] .> 1 ?
        # ijk[1] == 1 ? ione = 1 : ione = 0
        # ijk[2] == 1 ? jone = 1 : jone = 0
        # ijk[3] == 1 ? kone = 1 : kone = 0

        # +1 because of julia 1-indexing
        ione = (0, 1)[(p["gn"][1] > 1) + 1]
        jone = (0, 1)[(p["gn"][2] > 1) + 1]
        kone = (0, 1)[(p["gn"][3] > 1) + 1]


        if !p["guard_zones"] || !Bool(all)
            m = @SVector ones(Int32, 3)
            n = p["n"]
        elseif Bool(all)
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

    function indices_and_weights(p::Dict, point = [0.5, 0.5, 0.5], iv = 0)
        """
        Return indices and interpolation weights for a point in a patch p
        """

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


    function patches_along(pp, point=[0.5, 0.5, 0.5]; dir = 1, verbose=0)
        """
        Get the patches along a given direction through a point
        """

        "patches" in keys(pp) ? pp = pp["patches"] : nothing

        pt = copy(point)

        ## ADD VERBOSE ??? ###

        p = patch_at(pp, pt)

        if p != nothing
            p1 = p
            while p != nothing
                pt[dir] = p["position"][dir] - p["size"][dir]*0.6
                p = patch_at(pp, pt)
                p != nothing ? p1 = p : nothing
        end

        p = p1
        out = []
        while p != nothing
            push!(out, p)
            pt[dir] = p["position"][dir] + p["size"][dir]*0.6
            p = patch_at(pp, pt)
        end
        else
            out = nothing
        end

        return out

    end

    function patch_at(pp, point=[0.5, 0.5, 0.5]; verbose=0)
        """
        Find the patch that contains a given point
        """

        level = -1
        p1 = nothing
        @inbounds for p in pp
            if is_inside(p, point, verbose)
            p["level"] > level ? p1 = p : level = p["level"]
            end
        end

        return p1

    end

    function is_inside(p, point, verbose)
    """
    Return true/false depending on if the point is inside the patch
    """

    left = point .>= p["llc_cart"]
    right = point .<= p["llc_cart"] + p["size"]

    all(left) && all(right)

    end


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
    plane(patch::Dict; x::Float, y::Float, z::Float,
        iv::Union{String, Int}, all::Bool, verbose::Int)

    Return patch data of quantity iv at a slice x/y/z.

    Arguments:
    --------------
        - patch: Dictionairy, a patch object from a snapshot

    Kwargs:
    -------------
        - x, y, z: Float, position at which to slice, default nothing
        - iv: String/Int, what quantity to extract, default 0
        - all: Bool, whether to include guard zones, default false

    Returns:
    -------------
        - 2d array of Float32.
    """
    function plane(patch; x = nothing, y = nothing, z = nothing, iv = 0,
                        verbose = 0, all=false)



        if patch["guard_zones"]# && !all
            li = patch["li"]  # lower inner
            ui = patch["ui"]  # upper inner
        else #all || !patch["guard_zones"]
            li = ones(Int, 3) # using Static??
            ui = patch["n"]
        end

        #@assert !all(isnothing.((x, y, z))) "plane: must give one x, y, or z-value."

        variv = patch["var"](iv)
        if !isnothing(x)
            p = (x - patch["x"][1])/patch["ds"][1]
            i = round(Int, p)
            i = min(i, ui[1]-1)
            p -= i

            f = variv[i, li[2]:ui[2], li[3]:ui[3]]*(1.0 - p) +
                variv[i+1, li[2]:ui[2], li[3]:ui[3]]*p

        elseif !isnothing(y)
            p = (y - patch["y"][1])/patch["ds"][2]
            i = round(Int, p)
            i = min(i, ui[2]-1)
            p = p - i
            f = transpose(variv[li[1]:ui[1], i  , li[3]:ui[3]]*(1.0-p) +
                          variv[li[1]:ui[1], i+1, li[3]:ui[3]]*p)

        elseif !isnothing(z)
            p = (z - patch["z"][1])/patch["ds"][3]
            i = round(Int, p)
            i = min(i, ui[3]-1)
            p = p - i
            # if i == 0
            #     f = variv[li[1]:ui[1], li[2]:ui[2], end]*(1.0 - p) +
            #         variv[li[1]:ui[1], li[2]:ui[2], i+1]*p
            # else
            #     f = variv[li[1]:ui[1], li[2]:ui[2], i]*(1.0 - p) +
            #         variv[li[1]:ui[1], li[2]:ui[2], i+1]*p
            # end
            if i == 0
                f = variv[:, :, end]*(1.0 - p) +
                    variv[:, :, i+1]*p
            else
                f = variv[:, :, i]*(1.0 - p) +
                    variv[:, :, i+1]*p
            end
        end

        Bool(verbose) ? println("plane: $i, p = $i $p") : nothing

        return f

    end


    """
    box(patch::Dict; iv::Union{String, Int}, all::Bool, verbose::Int)

    Return volume of data of quantity iv from patch.

    Arguments:
    ------------
        - patch: Dictionairy, a patch object from a snapshot

    Kwargs:
    --------------
        - iv: String/Int, what quantity to extract, default 0
        - all: Bool, whether to include guard zones, default false

    Returns:
    ------------
        - 3d array of Float32.
    """
    function box(patch; iv = 0, all=false, verbose = 0)

        if patch["guard_zones"]# && !all
            li = patch["li"]  # lower inner
            ui = patch["ui"]  # upper inner
        else #all || !patch["guard_zones"]
            li = ones(Int, 3) # using Static??
            ui = patch["n"]
        end

        f = patch["var"](iv)#[li[1]:ui[1], li[2]:ui[2], li[3]:ui[3]]

        return f
    end


end #module