function image_plane(snap; x = nothing, y = nothing, z = 0.5, kw...)

     kw = Dict(kw)
     kv = Dict{Symbol, Any}(:verbose => 0, :iv => 0,
                            :grids => false, :cbar => false,
                            :title => 0, :label => nothing)


    _kw_extract(kw, kv)
    iv = kv[:iv]
    if kv[:label] == nothing
        kv[:label] = get_unit(snap, iv)
        if kv[:label] == nothing
            println("Unit for $iv not found")
            kv[:label] = string(iv)
        end
    end


    origin = snap["cartesian"]["origin"]
    sze = snap["cartesian"]["size"]

    center = origin + 0.5*sze

    i = nothing
    if x != nothing
        i = 1
        x = x*sze[i]
        center = center[[2, 3]]
        sze = sze[[2, 3]]
        dims = snap["cartesian"]["dims"][2:3]

    elseif y != nothing
        i = 3
        y = y*sze[i]
        center = center[[1, 3]]
        sze = sze[[1, 3]]
        dims = snap["cartesian"]["dims"][1:2:3]

    elseif z != nothing
        i = 3
        z = z*sze[i]
        center = center[[1, 2]]
        sze = sze[[1, 2]]
        dims = snap["cartesian"]["dims"][1:2]
    end

    xyz = (x, y, z)

    sizedir = ["x" "y" "z"]
    labels=[("y","z"), ("x","z"), ("x","y")]

    patches = patches_in(snap, x=x, y=y, z=z)

    ll = Dict{Int, Tuple{Array{Float64, 2}, Array{Float64, 1}, Int}}()

    p0 = patches[1]

    f = plane(p0, x=x, y=y, z=z, iv=iv, verbose=kv[:verbose])
    fmin = minimum(f)
    fmax = maximum(f)

    e = p0["extent"][i,:]

    emin = [e[1], e[3]]
    emax = [e[2], e[4]]

    xmin, xmax, ymin, ymax = (0, 0, 0, 0)
    datasize = 0
    for p in patches
        im = plane(p, x = x, y = y, z = z, iv = iv)
        datasize += prod(size(im))

        fmin = min(fmin, minimum(im))
        fmax = max(fmax, maximum(im))

        e = p["extent"][i,:]

        emin = min.(emin, [e[1], e[3]])
        emax = max.(emax, [e[2], e[4]])

        ll[p["id"]] = (im, e, p["level"])

        xmin = min(xmin, p["x"][1])
        xmax = max(xmax, p["x"][end])

        ymin = min(ymin, p["y"][1])
        ymax = max(ymax, p["y"][end])
    end

    fmin < 0 ? cmap = :coolwarm : cmap = :auto

    sortedkeys = sort(collect(keys(ll)))

    datashp = nothing
    try
        n = Int(sqrt(datasize))
        datashp = (n, n)
    catch
        datashp = getindex(p0["n"], [k for k in 1:3 if k != i]) .* dims
    end

    data = reshape(zeros(Float32, datasize), datashp...)

    for id in sortedkeys
        im, e = ll[id][1:2]

        xs = Int.(round.(e[1:2] .* datashp))
        ys = Int.(round.(e[3:4] .* datashp))

        data[(xs[1] + 1):xs[end], (ys[1] + 1):ys[end]] = im
    end

    xx = collect(LinRange(xmin, xmax, datashp[1]))
    yy = collect(LinRange(ymin, ymax, datashp[2]))

    hm = heatmap(xx, yy, data, colorbar_title=kv[:label], c=cmap)

    if kv[:title] != 0
        title!(kv[:title])
    else
        title!("$(sizedir[i]) = $(xyz[i])")
    end

    if kv[:grids]
        for p in patches
            e = p["extent"][i,:]

            x = [e[1], e[1], e[2], e[2], e[1]]
            y  =[e[3], e[4], e[4], e[3], e[3]]

            plot!(x, y, color=:white, label=false)
        end
    end

    xlabel!(labels[i][1])
    ylabel!(labels[i][2])

    return hm


end


function SlicePlot(snaps::Union{String, Tuple, AbstractArray}="all"
                       ; x = nothing, y = nothing, z = 0.5, kw...)

     kw = Dict(kw)
     kv = Dict{Symbol, Any}(:verbose => 0, :iv => 0,
                            :grids => false, :cbar => false,
                            :title => 0, :label => nothing,
                            :fps => 1, :outname => nothing,
                            :filetype => "mp4")


    _kw_extract(kw, kv)

    if kv[:outname] == nothing
        kv[:outname] = kv[:iv] * "." * kv[:filetype]
    end

    if typeof(snaps) <: String && snaps == "all"
        snapIDS = 0:get_n_snapshots()-1
        snaps = [snapshot(i) for i in snapIDS]
        snaps = snaps[snaps .!= nothing]
    else
        snaps = [snapshot(i) for i in snaps]
        snaps = snaps[snaps .!= nothing]
    end

    iv = kv[:iv]

    """ Plot first heatmap """
    snap0 = snaps[1]

    if kv[:label] == nothing
        kv[:label] = get_unit(snap0, iv)
        if kv[:label] == nothing
            println("Unit for $iv not found")
            kv[:label] = string(iv)
        end
    end


    origin = snap0["cartesian"]["origin"]
    sze = snap0["cartesian"]["size"]

    center = origin + 0.5*sze

    i = nothing
    if x != nothing
        i = 1
        x = x*sze[i]
        center = center[[2, 3]]
        sze = sze[[2, 3]]
        dims = snap0["cartesian"]["dims"][2:3]

    elseif y != nothing
        i = 3
        y = y*sze[i]
        center = center[[1, 3]]
        sze = sze[[1, 3]]
        dims = snap0["cartesian"]["dims"][1:2:3]

    elseif z != nothing
        i = 3
        z = z*sze[i]
        center = center[[1, 2]]
        sze = sze[[1, 2]]
        dims = snap0["cartesian"]["dims"][1:2]
    end

    xyz = (x, y, z)

    sizedir = ["x" "y" "z"]
    labels=[("y","z"), ("x","z"), ("x","y")]

    anim = @animate for snap in snaps
        ll = Dict{Int, Tuple{Array{Float64, 2}, Array{Float64, 1}, Int}}()
        patches = patches_in(snap, x=x, y=y, z=z)

        p0 = patches[1]

        f = plane(p0, x=x, y=y, z=z, iv=iv, verbose=kv[:verbose])
        fmin = minimum(f)
        fmax = maximum(f)

        e = p0["extent"][i,:]

        emin = [e[1], e[3]]
        emax = [e[2], e[4]]

        xmin, xmax, ymin, ymax = (0, 0, 0, 0)
        datasize = 0
        for p in patches
            im = plane(p, x = x, y = y, z = z, iv = iv)
            datasize += prod(size(im))

            fmin = min(fmin, minimum(im))
            fmax = max(fmax, maximum(im))

            e = p["extent"][i,:]

            emin = min.(emin, [e[1], e[3]])
            emax = max.(emax, [e[2], e[4]])

            ll[p["id"]] = (im, e, p["level"])

            xmin = min(xmin, p["x"][1])
            xmax = max(xmax, p["x"][end])

            ymin = min(ymin, p["y"][1])
            ymax = max(ymax, p["y"][end])
        end

        fmin < 0 ? cmap = :coolwarm : cmap = :auto

        sortedkeys = sort(collect(keys(ll)))

        data = reshape(zeros(Float32, datasize), Int(sqrt(datasize)), :)
        datashp = size(data)

        for id in sortedkeys
            im, e = ll[id][1:2]

            xs = Int.(round.(e[1:2] .* datashp))
            ys = Int.(round.(e[3:4] .* datashp))

            data[(xs[1] + 1):xs[end], (ys[1] + 1):ys[end]] = im
        end

        xx = collect(LinRange(xmin, xmax, datashp[1]))
        yy = collect(LinRange(ymin, ymax, datashp[2]))

        heatmap(xx, yy, data, colorbar_title=kv[:label], c=cmap)

        if kv[:grids]
            for p in patches
                e = p["extent"][i,:]

                x = [e[1], e[1], e[2], e[2], e[1]]
                y  =[e[3], e[4], e[4], e[3], e[3]]

                plot!(x, y, color=:white, label=false)
            end
        end

        if kv[:title] != 0
            title!(kv[:title] * " t = " * round(snap["time"], digits=2))
        else
            title!("t = $(round(snap["time"], digits=2))")
        end
        xlabel!(labels[i][1])
        ylabel!(labels[i][2])

    end
    gif(anim, kv[:outname], fps=kv[:fps]);

end

function test()
    z = 0.0039
    nsnaps = get_n_snapshots()

    snaps = [snapshot(i) for i = 0:nsnaps-1]
    snaps = snaps[snaps .!= nothing]

    data = [plane_buffer(snap, iv="b1", z=z) for snap in snaps]

    times = [snap["time"] for snap in snaps]

    t = range(times[1], times[end], length=2*length(times) - 1)
    newdata = [(data[i+1] .+ data[i])/2 for i = 1:length(data)-1]

    newwdata = []
    ctr1 = 1
    ctr2 = 1
    for k in 1:length(newdata)
        if (Bool(k%2))
            push!(newwdata, data[ctr1])
            ctr1 += 1
        else
            push!(newwdata, newdata[ctr2])
            ctr2 += 1
        end
    end


    anim = @animate for d in newwdata
        heatmap(d)
    end


    gif(anim, "out.mp4", fps=2)
end

function animate_surface(snaps::Union{String, Tuple, AbstractArray}="all"
                         ;x = nothing, y = nothing, z = 0.5, kw...)

    kw = Dict(kw)
    kv = Dict{Symbol, Any}(:verbose => 0, :iv => 0,
                           :cmap => nothing, :outpath => nothing,
                           :title => nothing, :label => nothing,
                           :fps => 2, :filetype => "mp4")

    _kw_extract(kw, kv)

    iv = kv[:iv]

    if typeof(snaps) <: String && snaps == "all"
        snapIDS = 0:get_n_snapshots()-1
        snaps = [snapshot(i) for i in snapIDS]
        snaps = snaps[snaps .!= nothing]
    else
        snaps = [snapshot(i) for i in snaps]
        snaps = snaps[snaps .!= nothing]
    end

    if kv[:label] == nothing
        kv[:label] = get_unit(snaps[1], iv)
    end

    xyz = [x, y, z]
    dirs = [(x, "x"), (y, "y"), (z, "z")]
    sliceDir = getindex(dirs, xyz .!= nothing)[1]
    planeDirs = getindex(dirs, xyz .== nothing)

    kv[:cmap] != nothing ? cmap = kv[:cmap] : cmap = :auto

    if typeof(kv[:title]) <: String
        title!(kv[:title])
    elseif kv[:title] == :pos
        title!("$(sliceDir[2]) = $(sliceDir[1])")
    elseif kv[:title] == :time
        time = round(snap["time"], digits=2)
        title!("t = $time")
    end


    anim = @animate for snap in snaps

        data = plane_buffer(snap, iv=iv, x=x, y=y, z=z)
        X = range(0, 1, length=size(data)[1])
        Y = range(0, 1, length=size(data)[2])
        surface(X, Y, data, label=kv[:label], c=cmap, cbar=false)

        xlabel!(planeDirs[1][2])
        ylabel!(planeDirs[2][2])

    end

    if kv[:outpath] != nothing
        outfile = kv[:outpath]*"."*kv[:filetype]
    else
        outfile = "surf_$iv."*kv[:filetype]
    end

    gif(anim, outfile, fps=kv[:fps])

    function Streamplot_(snap::Dict; iv = 0, x = nothing, y = nothing, z = 0.5,
                        unigrid=true, kw...)
        """ Plot a streamplot of variable iv over a heatmap/contour """

        kw = Dict(kw)
        kv = Dict{Symbol, Any}(:verbose => 0,
                               :grids => false, :cmap => nothing,
                               :title => nothing, :label => nothing,
                               :style => "heatmap", :fill => true,
                               :center => [0.5, 0.5], :width => 1.0,
                               :xlabel => nothing, :ylabel => nothing,
                               :savepath => nothing, :filetype => "png")

        _kw_extract(kw, kv)

        xyz = [x, y, z]
        dirs = [(x, "x"), (y, "y"), (z, "z")]
        sliceDir = getindex(dirs, xyz .!= nothing)[1]
        planeDirs = getindex(dirs, xyz .== nothing)

        fig = figure()
        ax = fig.gca()

        # quantity magnitude
        if unigrid
            d1, d2, Q = unigrid_plane(snap, iv=iv, x=x, y=y, z=z)
        else
            d1, d2, Q = amr_plane(snap, iv=iv, x=x, y=y, z=z)
        end

        kv[:transpose] ? Q = Q' : nothing

        # velocities
        u1 = zeros(Float32, size(Q))
        u2 = zeros(Float32, size(Q))

        u1[1:end-1,:] = diff(Q, dims=1)
        u2[:,1:end-1] = diff(Q, dims=2)

        nx, ny = size(Q)

        if kv[:style] == "heatmap"
            im = ax.imshow(Q, extent=[d1[1], d1[end], d2[1], d2[end]], origin="lower")
        elseif kv[:style] == "contour"
            im = ax.contourf(d1, d2, Q)
        end

        ax.streamplot(collect(d1), collect(d2), u1, u2, color="white", linewidth=0.2, density=1.5)


        if kv[:xlabel] != nothing
            ax.set_xlabel(kv[:xlabel])
        else
            ax.set_xlabel(planeDirs[1][2])
        end

        if kv[:ylabel] != nothing
            ax.set_ylabel(kv[:ylabel])
        else
            ax.set_ylabel(planeDirs[1][2])
        end

        if kv[:title] != nothing
            ax.set_title(kv[:title])
        end

        label = split(get_unit(snap, iv), "\$")
        quantLabel = label[1]
        unitLabel = latexstring(label[2])

        fig.colorbar(im, label=quantLabel * " magnitude " * unitLabel)

        if kv[:savepath] == nothing
            savepath = string(snap["iout"])*"_"*quantLabel*"_streamplot.png"
        else
            savepath = kv[:savepath]*"."*kv[:filetype]
        end
        #fig.savefig(savepath, dpi=300)

        return fig

    end
