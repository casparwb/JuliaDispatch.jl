using Plots, LaTeXStrings
using PyPlot: streamplot, imshow, figure
#import Makie
import WGLMakie
using JuliaDispatch.Utils, JuliaDispatch.Buffers
using JuliaDispatch.Select
using JuliaDispatch.Dispatch: snapshot
using Unitful, Latexify, UnitfulRecipes, UnitfulLatexify
WGLMakie.activate!()
gr()

default(:size, (1200, 800))


""" 
    plot_values_along(snap::Dict, pt::Array=[0.5, 0.5, 0.5]; iv=0, dir=1, kw...)

Plot values of quantity `iv` along direction `dir={1,2,3}` (x, y, z), through point `pt=[x,y,z]`. 

Accepts any keyword arguments supported by the `plot`-function exported from Plots.GR(). See 
[documentation](https://docs.juliaplots.org/latest/generated/gr/) for a complete list.

# Examples
``` 
julia> plot_values_along(snap, iv="ekin", linestyle=:scatter, xlabel="x", ylabel="y")
julia> plot_values_along(snap, iv="T", dir=3, style=:dot, label="Temperature")
julia> plot_values_along(snap, [0.0, 1.0, 2.0], iv="ekin", title="Kinetic Energy")
```
"""
function plot_values_along(snap::Dict, pt=[0.5, 0.5, 0.5]; kw...)
    kw = Dict(kw)
    # kv = Dict{Symbol, Any}(:dir => 1, :verbose => 0,
    #                        :all => false, :iv => 0,
    #                        :i4 => 0, :label => false,
    #                        :ylabel => 0, :xlabel => 0,
    #                        :title => 0, :style=> "line",
    #                        :norm => false, :bins => :auto)

    kv = Dict{Symbol, Any}(:verbose => 0,
                           :all => false, :iv => 0,
                           :i4 => 0)

    _kw_extract(kw, kv)

    data = values_along(snap, pt, iv=kv[:iv], dir=kv[:dir], all=kv[:all], verbose=kv[:verbose])

    # if kv[:ylabel] == 0
    #     kv[:ylabel] = get_unit(snap, kv[:iv])
    # end

    time = round(snap["time"], digits=2)
    # plt = plot(data, label=kv[:label], st=kv[:style],
    #            normalize=kv[:norm], bins=kv[:bins])
    plt = plot(data; kw...)
    # if kv[:xlabel] == 0
    #     if kv[:dir] == 1
    #         dir = "x"
    #         xlabel!("x")
    #     elseif kv[:dir] == 2
    #         dir = "y"
    #         xlabel!("y")
    #     elseif kv[:dir] == 3
    #         dir = "z"
    #         xlabel!("z")
    #     end
    # else
    #     xlabel!(kv[:xlabel])
    # end

    # # ylabel!(kv[:ylabel])

    # if kv[:title] != 0
    #     title!(kv[:title])
    # end

    # plt

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


"""
    sliceplot(snap::Dict; x = nothing, y = nothing, z = nothing, unigrid = true, kw...)


Make a sliceplot of quantity ´iv´ at the given coordinate. Accepts any keyword argument
supported by the `Plots.plot()` function, in addition to

# Arguments
- `iv::Union{Int, String, Char}`: what quantity to plot. Default `0`.
- `grids::Bool`: whether to show patch grids. Default `false`.
- `width::Tuple{Int, Int}`: width of the sub-domain axes. Default `nothing` (whole slice is plotted). 
- `center::Tuple{Int, Int}`:  where to center the plot. Default `nothing` (center is not moved)
- `resample::Union{Int, Tuple}`: whether to resample (upscale or downscale) slice. Default `nothing`. Can
                                  be `Int` or `Tuple` of new dimensions.
- `dims::Union{Tuple{Int, Int}, Int}`: size of interpolated slice if experiment uses adaptive mesh-refinement. 
                                       Can be either tuple indicating size in each dimension, or
                                       an integer indicating same size for both dimensions. Default `300`.


# Examples
```
```
"""
function sliceplot(snap::Dict,
                  ;x = nothing, y = nothing, z = nothing, unigrid=true,
                  kw...)

    kw = Dict{Symbol, Any}(kw)
    if !haskey(kw, :linetype) 
        kw[:linetype] = :heatmap 
    end
    kv = Dict{Symbol, Any}(:verbose => 0, :iv => 0,
                           :grids => false,
                           :width => nothing, :dims => 300,
                           :center => nothing, :resample => nothing,
                           :resampledims => nothing)

    _kw_extract(kw, kv)
    iv = kv[:iv]
    verbose = kv[:verbose]

    # get the axis where to slice and the resulting plane
    xyz = [x, y, z]
    dirs = [(x, "x", 1), (y, "y", 2), (z, "z", 3)]
    axis = getindex(dirs, xyz .!= nothing)[1]
    planeDirs = getindex(dirs, xyz .== nothing)

    verbose > 1 && println("Sliceplot through $(axis[2]) in plane $((planeDirs[1][2], planeDirs[2][2]))")

    origin = copy(snap["cartesian"]["origin"])
    Size = copy(snap["cartesian"]["size"])
    deleteat!(origin, axis[3])
    deleteat!(Size, axis[3])

    if unigrid
        data = unigrid_plane(snap, iv=iv, x=x, y=y, z=z)
        verbose >= 1 && println("Unigrid data with shape $(size(data))")
        println(dimension(eltype(data)))

    else
        data = amr_plane(snap, iv=iv, x=x, y=y, z=z, dims=kv[:dims])
        data = data'
        verbose >= 1 && println("Mesh refined data with shape $(size(data))")
    end

    # plane axes
    d1 = range(origin[1], origin[1]+Size[1], length=size(data)[1])
    d2 = range(origin[2], origin[2]+Size[2], length=size(data)[2])

    # check if width is given
    wflag = false
    if kv[:width] !== nothing
        width = kv[:width]
        wflag = true
    else
        width = Size
    end

    # check if center is given
    cflag = false
    if kv[:center] !== nothing
        center = kv[:center]
        cflag = true
    else
        center = origin .+ width
    end

    # realign data if new width of center is given
    if wflag || cflag
        idxs1 = findall(abs.(d1 .- center[1]) .<= width[1])
        idxs2 = findall(abs.(d2 .- center[2]) .<= width[2])

        data = data[idxs1, idxs2]
        d1 = d1[idxs1]
        d2 = d2[idxs2]

    end

    # resample data if given
    if !isnothing(kv[:resample])
        newdims =  kv[:resample]
        verbose >= 1 && println("Resampling to size $newdims")
        d1, d2, data = resample(d1, d2, data, newdims)
    end

    #hm = plot(d1, d2, data'; kw...)
    hm = plot(data; unitformat=latexify, kw...)
    if kv[:grids]
        i = axis[3]
        patches = patches_in(snap, x=x, y=y, z=z)
        for p in patches
            e = p["extent"][i,:]

            x = [e[1], e[1], e[2], e[2], e[1]]
            y = [e[3], e[4], e[4], e[3], e[3]]

            plot!(hm, x, y, color=:gray, label=false)
        end
    end

    # if kv[:title] == ""
    #     pos = (axis[2], axis[1])
    #     time = round(snap["time"], digits=2)
    #     title!("$(pos[1]) = $(pos[2]), t = $time, iv=$iv")
    # elseif typeof(kv[:title]) <: String
    #     title!(kv[:title])
    # end


    # if kv[:xlabel] !== nothing
    #     xlabel = kv[:xlabel]
    # else
    #     xlabel = planeDirs[1][2]
    # end

    # if kv[:ylabel] !== nothing
    #     ylabel = kv[:ylabel]
    # else
    #     ylabel = planeDirs[2][2]
    # end


    # xlabel!(xlabel)
    # ylabel!(ylabel)

    return hm

end

function test_heatmap(data)
    heatmap(data, clims=(0.0u"m/s", 2.0u"m/s"))
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

function volume(snap::Dict; iv = 0, unigrid=true, kw...)
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

    if unigrid
        data = unigrid_volume(snap, iv=iv)
    else
        data = amr_volume(snap, iv=iv)
    end
    
    nx, ny, nz = size(data)

    x = range(start[1], stop[1], length=nx)
    y = range(start[2], stop[2], length=ny)
    z = range(start[3], stop[3], length=nz)


    # kv[:cmap] !== nothing ? cmap = kv[:cmap] : cmap = :viridis
    scene = WGLMakie.volume(x, y, z, data)
    # kv[:title] !== nothing ? Makie.title(scene, kv[:title]) : nothing

    # if kv[:grids]
    #     for patch in snap["patches"]
    #         e = patch["extent"]
    #     end
    # end
    WGLMakie.save("volume.png", scene)
    #scene
end

function anim_plane(;data="../data", run="", x = nothing, y = nothing, z = nothing, iv=0, 
                    tspan=nothing, unigrid=true, step=1, savepath=nothing)

    nsnapshots = get_n_snapshots(data=data, run=run)

    anim = @animate for i = 1:step:nsnapshots-1
        snap = snapshot(i, data=data, run=run)
        isnothing(snap) && continue
        plne = unigrid_plane(snap, x=x, y=y, z=z, iv=iv)
        heatmap(plne)
    end #every step

    gif(anim, savepath)
    

end

function animate_volume(;data="../data", run="", iv=0, unigrid=true, savepath=".")

    if unigrid
        data = unigrid_volume(snap, iv=iv)
    else
        data = amr_volume(snap, iv=iv)
    end


    nsnaps = get_n_snapshots(data=data, run=run)
    snap0 = snapshot(0, data=data, run=data)

    scene = WGLMakie.volume(data)
    vol_plot = scene[end]

    record(savepath, 2:nsnaps-1) do i
        vol_plot[1] = unigrid_volume(snapshot(i, data=data, run=run), iv=iv)
    end

end