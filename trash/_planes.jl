function amr_plane(snap; iv = 0, x = nothing, y = nothing, z = 0.5, Log = false,
                   grids=false)

    xyz = [x, y, z]
    dir = getindex((1, 2, 3), xyz .!= nothing)[1]
    ax = xyz[xyz .!= nothing][1]

    patches = patches_in(snap, x=x, y=y, z=z)

    fig = figure(figsize=(12, 16))
    ax = fig.gca()

    emin = [Inf, Inf]
    emax = [-Inf, -Inf]

    vmin = Inf
    vmax = -Inf

    patchDict = Dict()

    for patch in patches
        e = patch["extent"][dir,:]
        data = interpolate(patch, iv=iv, x=x, y=y, z=z)

        maxValue = maximum(data)
        minValue = minimum(data)

        if maxValue > vmax
            vmax = maxValue
        end

        if minValue < vmin
            vmin = minValue
        end

        emin[1] = min(emin[1], e[1])
        emin[2] = min(emin[2], e[3])

        emax[1] = max(emax[1], e[2])
        emax[2] = max(emax[2], e[4])

        patchDict[patch["id"]] = (e, data)

    end

    e, data = patchDict[1]
    im = ax.imshow(data', extent=e, origin="lower", interpolation="nearest",
                vmin=vmin, vmax=vmax)

    for id in sort(collect(keys(patchDict)))[2:end]
        e, data = patchDict[id]
        Log ? data = log10.(data) : nothing

        ax.imshow(data', extent=e, origin="lower", interpolation="nearest",
                       vmin=vmin, vmax=vmax)

        if grids
            # d1, d2 = outline(e)
            d1 = [e[1], e[1], e[2], e[2], e[1]]
            d2 = [e[3], e[4], e[4], e[3], e[3]]
            ax.plot(d1, d2, color="white")
        end

        # ax.text(e[1], e[3], "$id")
    end

    ax.set_xlim([emin[1], emax[1]])
    ax.set_ylim([emin[2], emax[2]])


    return fig, im

end

function amr_volume(snap; iv::Union{Int, Array, String} = 0, all = true,
                    dims::Union{Tuple, Int}=100, verbose=0)
    """
    function for creating a unigrid volume buffer of 3d amr data
    """

    if typeof(dims) <: Tuple
        nx, ny, nz = dims
    else
        nx = ny = nz = dims
    end

    buffer = nothing
    ivs = nothing
    if typeof(iv) <: String && iv == "all"
        ivs = [iv for (k, iv) in snap["idx"] if (typeof(k) <: String && iv > 0)]
        buffer = Dict{Union{String, Int}, Array{Float32, 3}}()
        for iv in ivs
            buffer[iv] = zeros(Float32, nx, ny, nz)
        end
    elseif typeof(iv) <: Array
        ivs = iv
        buffer = Dict{Union{String, Int}, Array{Float32, 3}}()
        for iv_ in iv
            buffer[iv_] = zeros(Float32, nx, ny, nz)
        end
    else
        ivs = [iv]
        buffer = Dict{Union{String, Int}, Array{Float32, 3}}
        buffer[iv] = zeros(Float32, nx, ny, nz)
    end

    patches = snap["patches"]

    Size = snap["cartesian"]["size"]
    origin = snap["cartesian"]["origin"]
    periodic = snap["periodic"]

    x = collect(range(origin[1], origin[1]+Size[1], length=nx))
    y = collect(range(origin[2], origin[2]+Size[2], length=ny))
    z = collect(range(origin[3], origin[3]+Size[3], length=nz))
    for iv in ivs
        for patch in patches
            data = box(patch, iv=iv, verbose=verbose, with_gz=true)

            e = patch["extent"]
            ex = e[3,1:2]
            ey = e[3,3:4]
            ez = e[2,1:2]

            xs = findall(ex[1] .<= x .<= ex[2])
            ys = findall(ey[1] .<= y .<= ey[2])
            zs = findall(ez[1] .<= z .<= ez[2])

            @inbounds for ix in xs, iy in ys, iz in zs
                buffer[iv][ix, iy, iz] = interpolate_3d(patch, data,
                                                    x = x[ix],
                                                    y = y[iy],
                                                    z = z[iz],
                                                    periodic=periodic)

            end

        end
    end

    if length(keys(buffer)) == 1
        return collect(values(buffer))[1]
    else
        return buffer
    end
end
