using Plots, LaTeXStrings
using PyPlot: streamplot, imshow, figure
import Makie
using JuliaDispatch.Utils, JuliaDispatch.Buffers
using JuliaDispatch.Select
using JuliaDispatch.Dispatch: snapshot
gr()
default(:size, (1200, 800))


#include("../select/_select.jl")
#include("utils.jl")
#include("../select/buffers.jl")


function plot_values_along(snap::Dict, pt=[0.5, 0.5, 0.5]; kw...)
    """ Plot values along direction dir={1,2,3}, through point pt=[x,y,z] """
    kw = Dict(kw)
    kv = Dict{Symbol, Any}(:dir => 1, :verbose => 0,
                           :all => false, :iv => 0,
                           :i4 => 0, :label => false,
                           :ylabel => 0, :xlabel => 0,
                           :title => 0, :style=> "line",
                           :norm => false, :bins => :auto)

    _kw_extract(kw, kv)

    data = values_along(snap, pt, iv=kv[:iv], dir=kv[:dir], all=kv[:all])

    if kv[:ylabel] == 0
        kv[:ylabel] = get_unit(snap, kv[:iv])
    end

    time = round(snap["time"], digits=2)
    plt = plot(data, label=kv[:label], st=kv[:style],
               normalize=kv[:norm], bins=kv[:bins])

    if kv[:xlabel] == 0
        if kv[:dir] == 1
            dir = "x"
            xlabel!("x")
        elseif kv[:dir] == 2
            dir = "y"
            xlabel!("y")
        elseif kv[:dir] == 3
            dir = "z"
            xlabel!("z")
        end
    else
        xlabel!(kv[:xlabel])
    end

    # ylabel!(kv[:ylabel])

    if kv[:title] != 0
        title!(kv[:title])
    end

    plt

end


function plot_density(snap; iv=0)

end

function plot_stacked_density(snap; iv=0)

end

function histogram_along(snap, pt=[0.5, 0.5, 0.5]; kw...)
    """ Plot a histogram along direction dir={1,2,3}, through point pt=[x,y,z]
    kwargs
    :


    """

    kw = Dict(kw)
    kv = Dict{Symbol, Any}(:dir => 1, :verbose => 0,
                                                :all => false, :iv => 0,
                                                :i4 => 0, :var => 0,
                                                :label => nothing, :title => 0,
                                                :norm => false, :bins => :auto)


    _kw_extract(kw, kv)

    data = values_along(snap, pt, iv=kv[:iv], dir=kv[:dir], all=kv[:all])

    if kv[:label] == nothing
        kv[:label] = get_unit(snap, kv[:iv])
    end

    xyz = ["x", "y", "z"]
    idxs = [k for k in 1:3 if k != kv[:dir]]
    at = getindex(pt, idxs)
    at_idx = getindex(xyz, idxs)
    time = round(snap["time"], digits=2)

    label = "t = $time at $(at_idx[1]) = $(at[1]), $(at_idx[2]) = $(at[2])"

    hist = histogram(data, normalize=kv[:norm], label=label, bins=kv[:bins])
    if kv[:norm] == false
        ylabel!("count")
    else
        key = kv[:norm]
        label = (key == :pdf ? "probability density" : string(key))
        println(label)
        ylabel!(label)
    end


    xlabel!(kv[:label])

    if kv[:title] != 0
        title!(kv[:title])
    else
        title!("direction : $(xyz[kv[:dir]])")
    end

    hist
end
#
# function histogram_along!(snap, pt=[0.5, 0.5, 0.5]; kw...)
#     """
#     append to a histogram
#     """
#
#     kw = Dict(kw)
#     kv = Dict{Symbol, Any}(:dir => 1, :verbose => 0,
#                            :all => false, :iv => 0,
#                            :i4 => 0, :var => 0,
#                            :label => nothing, :title => 0,
#                            :norm => false, :bins => :auto)
#
#
#     _kw_extract(kw, kv)
#
#     data = values_along(snap, pt, iv=kv[:iv], dir=kv[:dir], all=kv[:all])
#     xyz = ["x", "y", "z"]
#     idxs = [k for k in 1:3 if k != kv[:dir]]
#     at = getindex(pt, idxs)
#     at_idx = getindex(xyz, idxs)
#     time = round(snap["time"], digits=2)
#
#     label = "t = $time at $(at_idx[1]) = $(at[1]), $(at_idx[2]) = $(at[2])"
#     histogram!(data, normalize=kv[:norm], label=label)
# end



# function Surface(snap::Dict; x = nothing, y = nothing, z = 0.5, unigrid=true,
#                 kw...)
#
#     kw = Dict(kw)
#     kv = Dict{Symbol, Any}(:verbose => 0, :iv => 0,
#                            :grids => false, :cmap => nothing,
#                            :title => nothing, :label => nothing,
#                            :xlims => (0, 1), :ylims => (0, 1),
#                            :xlabel => nothing, :ylabel => nothing)
#
#     _kw_extract(kw, kv)
#
#     iv = kv[:iv]
#
#     if kv[:label] == nothing
#         kv[:label] = get_unit(snap, iv)
#     end
#
#     xyz = [x, y, z]
#     dirs = [(x, "x"), (y, "y"), (z, "z")]
#     sliceDir = getindex(dirs, xyz .!= nothing)[1]
#     planeDirs = getindex(dirs, xyz .== nothing)
#
#     if unigrid
#         data = unigrid_plane(snap, x=x, y=y, z=z, iv=iv, verbose=verbose)
#     else
#         data = amr_plane(snap, x=x, y=y, z=z, iv=iv, verbose=verbose)
#     end
#
#     x0, xend = kv[:xlims]
#     y0, yend = kv[:ylims]
#
#     X = range(x0, xend, length=size(data)[1])
#     Y = range(y0, yend, length=size(data)[2])
#
#     kv[:cmap] != nothing ? cmap = kv[:cmap] : cmap = :auto
#
#     surface(X, Y, data, label=kv[:label], c=cmap, cbar=false)
#
#     if typeof(kv[:title]) <: String
#         title!(kv[:title])
#     elseif kv[:title] == :pos
#         title!("$(sliceDir[2]) = $(sliceDir[1])")
#     elseif kv[:title] == :time
#         time = round(snap["time"], digits=2)
#         title!("t = $time")
#     end
#
#     if kv[:xlabel] != nothing
#         xlabel!(kv[:xlabel])
#     else
#         xlabel!(planeDirs[1][2])
#     end
#
#     if kv[:ylabel] != nothing
#         ylabel!(kv[:ylabel])
#     else
#         ylabel!(planeDirs[2][2])
#     end
#
#
# end

function sliceplot(d1, d2, data::Array{T, 2} where T;
                  x = nothing, y = nothing, z = nothing, kw...)
    """
    Sliceplot wrapper for arrays
    """

    kw = Dict(kw)
    kv = Dict{Symbol, Any}(:verbose => 0,
                           :grids => false, :cmap => :auto,
                           :title => "", :label => nothing,
                           :style => :heatmap, :fill => true,
                           :width => nothing, :dims => 300,
                           :xlabel => nothing, :ylabel => nothing,
                           :center => nothing, :resample => false,
                           :resampledims => nothing)

    _kw_extract(kw, kv)
    verbose = kv[:verbose]

    wflag = false
    if kv[:width] != nothing
        width = kv[:width]
        wflag = true
    else
        width = [d1[end] - d1[1], d2[end] - d1[1]]
    end

    cflag = false
    if kv[:center] != nothing
        center = kv[:center]
        cflag = true
    else
        cidx1 = Int(round(length(d1)/2))
        cidx2 = Int(round(length(d2)/2))
        center = [d1[cidx1], d2[cidx2]]
    end

    if wflag || cflag
        idxs1 = findall(abs.(d1 .- center[1]) .<= width[1])
        idxs2 = findall(abs.(d2 .- center[2]) .<= width[2])

        data = data[idxs1, idxs2]
        d1 = d1[idxs1]
        d2 = d2[idxs2]
    end

    if kv[:resample]
        newdims = nothing
        if kv[:resampledims] == nothing
            newdims = kv[:dims]
        else
            newdims = kv[:resampledims]
        end
        Bool(verbose) ? println("Resampling with size $newdims") : nothing
        d1, d2, data = resample(d1, d2, data, newdims)
    end

    if kv[:label] == nothing
        kv[:label] = ""
    end

    hm = nothing
    if string(kv[:style]) == "streamplot"
        return Streamplot(d1, d2, data)
    else
        hm = plot(d1, d2, data,
                  st=kv[:style],
                  colorbar_title=kv[:label],
                  c=kv[:cmap],
                  fill=kv[:fill])
    end

    # if kv[:grids]
    #     i = axis[3]
    #     patches = patches_in(snap, x=x, y=y, z=z)
    #     for p in patches
    #         e = p["extent"][i,:]
    #
    #         x = [e[1], e[1], e[2], e[2], e[1]]
    #         y = [e[3], e[4], e[4], e[3], e[3]]
    #
    #         plot!(x, y, color=:white, label=false)
    #     end
    # end

    if kv[:title] != ""
        title!(kv[:title])
    end

    if kv[:xlabel] != nothing
        xlabel = kv[:xlabel]
    else
        xlabel = ""
    end

    if kv[:ylabel] != nothing
        ylabel = kv[:ylabel]
    else
        ylabel = ""
    end

    xlabel!(xlabel)
    ylabel!(ylabel)

    return hm

end

function sliceplot(snap::Dict,
                  ;x = nothing, y = nothing, z = nothing, unigrid=true,
                  kw...)
    """
    Make a sliceplot at the given coordinate.
    Input:
    -----------------
    snap: dict, snapshot
    x, y, z: float, at what coordinate to slice
    """

    kw = Dict(kw)
    kv = Dict{Symbol, Any}(:verbose => 0, :iv => 0,
                           :grids => false, :cmap => :auto,
                           :title => "", :label => nothing,
                           :style => :heatmap, :fill => true,
                           :width => nothing, :dims => 300,
                           :xlabel => nothing, :ylabel => nothing,
                           :center => nothing, :resample => false,
                           :resampledims => nothing)

    _kw_extract(kw, kv)
    iv = kv[:iv]
    verbose = kv[:verbose]

    xyz = [x, y, z]
    dirs = [(x, "x", 1), (y, "y", 2), (z, "z", 3)]
    axis = getindex(dirs, xyz .!= nothing)[1]
    planeDirs = getindex(dirs, xyz .== nothing)

    origin = copy(snap["cartesian"]["origin"])
    Size = copy(snap["cartesian"]["size"])
    deleteat!(origin, axis[3])
    deleteat!(Size, axis[3])

    if unigrid
        data = unigrid_plane(snap, iv=iv, x=x, y=y, z=z)
        Bool(verbose) ? println("Unigrid data with shape $(size(data))") : nothing
    else
        data = amr_plane(snap, iv=iv, x=x, y=y, z=z, dims=kv[:dims])
        data = data'
        Bool(verbose) ? println("Mesh refined data with shape $(size(data))") : nothing
    end

    d1 = range(origin[1], origin[1]+Size[1], length=size(data)[1])
    d2 = range(origin[2], origin[2]+Size[2], length=size(data)[2])

    wflag = false
    if kv[:width] != nothing
        width = kv[:width]
        wflag = true
    else
        width = Size
    end

    cflag = false
    if kv[:center] != nothing
        center = kv[:center]
        cflag = true
    else
        center = origin .+ width
    end

    if wflag || cflag
        idxs1 = findall(abs.(d1 .- center[1]) .<= width[1])
        idxs2 = findall(abs.(d2 .- center[2]) .<= width[2])

        data = data[idxs1, idxs2]
        d1 = d1[idxs1]
        d2 = d2[idxs2]

    end

    if kv[:resample]
        newdims = nothing
        if kv[:resampledims] == nothing
            newdims = kv[:dims]
        else
            newdims = kv[:resampledims]
        end
        Bool(verbose) ? println("Resampling with size $newdims") : nothing
        d1, d2, data = resample(d1, d2, data, newdims)
    end

    # if isnothing(kv[:label], nothing)
    #     kv[:label] = get_unit(snap, iv)
    # end

    hm = nothing
    if string(kv[:style]) == "streamplot"
        return Streamplot(d1, d2, data)
    else
        hm = plot(d1, d2, data,
                  st=kv[:style],
                  colorbar_title=kv[:label],
                  c=kv[:cmap],
                  fill=kv[:fill])
    end

    if kv[:grids]
        i = axis[3]
        patches = patches_in(snap, x=x, y=y, z=z)
        for p in patches
            e = p["extent"][i,:]

            x = [e[1], e[1], e[2], e[2], e[1]]
            y = [e[3], e[4], e[4], e[3], e[3]]

            plot!(x, y, color=:gray, label=false)
        end
    end

    if kv[:title] == ""
        pos = (axis[2], axis[1])
        time = round(snap["time"], digits=2)
        title!("$(pos[1]) = $(pos[2]), t = $time, iv=$iv")
    elseif typeof(kv[:title]) <: String
        title!(kv[:title])
    end


    if kv[:xlabel] !== nothing
        xlabel = kv[:xlabel]
    else
        xlabel = planeDirs[1][2]
    end

    if kv[:ylabel] !== nothing
        ylabel = kv[:ylabel]
    else
        ylabel = planeDirs[2][2]
    end


    xlabel!(xlabel)
    ylabel!(ylabel)

    return hm

end


function streamplot_(d1, d2, data)
    fig = figure(figsize=(12, 8))
    ax = fig.gca()

    u1 = zeros(Float32, size(data))
    u2 = zeros(Float32, size(data))

    u1[1:end-1,:] = diff(data, dims=1)
    u2[:,1:end-1] = diff(data, dims=2)

    # np = pyimport("numpy")
    # D1, D2 = np.meshgrid(d1, d2)

    im = ax.imshow(data, extent=[d1[1], d1[end], d2[1], d2[end]], origin="lower")
    ax.streamplot(collect(d1), collect(d2), u1, u2, color="white", linewidth=0.2, density=1.5)
    fig.colorbar(im)
    return fig
end

function anim_pane(snap; ax=1, nframes=10, unigrid=true, kw...)
    """ Animate a pane through given axis """

    kw = Dict(kw)
    kv = Dict{Symbol, Any}(:verbose => 0, :iv => 0,
                           :grids => false, :cmap => :auto,
                           :title => "", :label => nothing,
                           :style => :heatmap, :fill => true,
                           :width => nothing, :dims => 300,
                           :xlabel => nothing, :ylabel => nothing,
                           :center => nothing, :resample => false,
                           :resampledims => nothing)

    _kw_extract(kw, kv)
    iv = kv[:iv]
    verbose = kv[:verbose]

    axStr = ("x", "y", "z")[ax .== [1, 2, 3]][1]
    # planeDirs = getindex(dirs, xyz .== nothing)

    origin = copy(snap["cartesian"]["origin"])
    Size = copy(snap["cartesian"]["size"])
    deleteat!(origin, ax)
    deleteat!(Size, ax)

    get_data() = nothing
    # if unigrid
    #     if ax == 1
    #         get_data(xx) = unigrid_plane(snap, iv=iv, x=xx)
    #     elseif ax == 2
    #         get_data(yy) = unigrid_plane(snap, iv=iv, y=yy)
    #     else
    #         get_data(zz) = unigrid_plane(snap, iv=iv, z=zz)
    #     end
    # else
    #     if ax == 1
    #         get_data(xx) = amr_plane(snap, iv=iv, x=xx, y=nothing, z=nothing, dims=kv[:dims])'
    #     elseif ax == 2
    #         get_data(yy) = amr_plane(snap, iv=iv, y=yy,x=nothing, z=nothing, dims=kv[:dims])'
    #     else
    #         get_data(zz) = amr_plane(snap, iv=iv, z=zz,x=nothing, y=nothing, dims=kv[:dims])'
    #     end
    # end

    get_data(xx) = amr_plane(snap, iv=iv, x=xx, y=nothing, z=nothing, dims=kv[:dims])'

    start = snap["cartesian"]["origin"][ax]
    stop = snap["cartesian"]["size"][ax]

    axvals = range(start, start+stop, length=nframes+1)
    anim = @animate for i = 1:nframes
        data = get_data(axvals[i])
        d1 = range(origin[1], origin[1]+Size[1], length=size(data)[1])
        d2 = range(origin[2], origin[2]+Size[2], length=size(data)[2])

        wflag = false
        if kv[:width] != nothing
            width = kv[:width]
            wflag = true
        else
            width = Size
        end

        cflag = false
        if kv[:center] != nothing
            center = kv[:center]
            cflag = true
        else
            center = origin .+ width
        end

        if wflag || cflag
            idxs1 = findall(abs.(d1 .- center[1]) .<= width[1])
            idxs2 = findall(abs.(d2 .- center[2]) .<= width[2])

            data = data[idxs1, idxs2]
            d1 = d1[idxs1]
            d2 = d2[idxs2]

        end

        if kv[:resample]
            newdims = nothing
            if kv[:resampledims] == nothing
                newdims = kv[:dims]
            else
                newdims = kv[:resampledims]
            end
            Bool(verbose) ? println("Resampling with size $newdims") : nothing
            d1, d2, data = resample(d1, d2, data, newdims)
        end

        heatmap(d1, d2, data, cbar=false)
    end

    # gif(anim, "testgif.gif", fps=1)
    gif(anim, "testgif2.gif", fps=20)
end

function volume(snap::Dict; iv = 0, kw...)
    """ Plot a 3D volume of the given quantity iv """

    kw = Dict(kw)
    kv = Dict{Symbol, Any}(:verbose => 0, :iv => 0,
                           :grids => false, :cmap => nothing,
                           :title => nothing, :label => nothing,
                           :style => :heatmap, :fill => true,
                           :width => (1.0, 1.0),
                           :xlabel => nothing, :ylabel => nothing,
                           :transpose => false)

    _kw_extract(kw, kv)

    origin = snap["cartesian"]["origin"]
    size_ = snap["cartesian"]["size"]

    start = origin .- size_/2
    stop = origin .+ size_/2

    if typeof(iv) <: AbstractArray
        data = iv
    else
        data = box_buffer(snap, iv=iv)
    end

    nx, ny, nz = size(data)

    x = range(start[1], stop[1], length=nx)
    y = range(start[2], stop[2], length=ny)
    z = range(start[3], stop[3], length=nz)


    kv[:cmap] !== nothing ? cmap = kv[:cmap] : cmap = :viridis
    scene = Makie.volume(x, y, z, data, colormap=cmap, isovalue=1)
    kv[:title] !== nothing ? Makie.title(scene, kv[:title]) : nothing

    # if kv[:grids]
    #     for patch in snap["patches"]
    #         e = patch["extent"]
    #     end
    # end

    scene
end

function anim_plane(;data="", x = nothing, y = nothing, z = nothing, iv=0, 
                    tspan=nothing, unigrid=true, step=1)

    nsnapshots = get_n_snapshots(data=data)

    anim = @animate for i = 1:nsnapshots-1
        snap = snapshot(i, data=data)
        isnothing(snap) && continue
        plne = unigrid_plane(snap, x=x, y=y, z=z, iv=iv)
        heatmap(plne)
    end every step

    gif(anim, "test.gif")
    

end

