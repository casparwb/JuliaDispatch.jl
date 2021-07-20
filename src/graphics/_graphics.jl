using Plots, LaTeXStrings, ProgressBars
# using PyPlot: streamplot, imshow, figure
#import Makie
import WGLMakie
using JuliaDispatch.Utils, JuliaDispatch.Buffers
using JuliaDispatch.Select
using JuliaDispatch.Dispatch: snapshot
using Latexify, LaTeXStrings
# using Unitful, Latexify, UnitfulRecipes, UnitfulLatexify
WGLMakie.activate!()
gr()
default(:size, (800, 600))


""" 
    plot_values_along(snap::Dict, pt::Array=[0.5, 0.5, 0.5]; iv = 0, dir = 1, verbose = 0, kw...)

Plot values of quantity `iv` along direction `dir={1,2,3}` (x, y, z), through point `pt=[x,y,z]`. 

Accepts any keyword arguments supported by the `plot`-function exported from Plots.GR(). See 
[documentation](https://docs.juliaplots.org/latest/generated/gr/) for a complete list.

# Examples
``` 
julia> plot_values_along(snap, iv="ekin", linetype=:scatter, xlabel="x", ylabel="y") # scatter plot
julia> plot_values_along(snap, iv="T", dir=3, style=:dashdot, label="Temperature")
julia> plot_values_along(snap, [0.0, 1.0, 2.0], iv="ekin", title="Kinetic Energy")
```
"""
function plot_values_along(snap::Dict, pt=[0.5, 0.5, 0.5]; iv = 0, dir = 1, verbose = 0, kw...)
    kw = Dict(kw)


    kv = Dict{Symbol, Any}(:verbose => 0,
                           :all => false, :iv => 0,
                           :i4 => 0, :label => nothing)

    _kw_extract(kw, kv)
    pt = Float64.(pt)
    data = values_along(snap, pt, iv=iv, dir=dir, all=kv[:all], verbose=verbose)


    plt = plot(data, label=kv[:label]; kw...)

    if !haskey(kw, :xlabel)
        xlabel!(plt, [L"x", L"y", L"z"][dir])
    end

    if !haskey(kw, :ylabel)
        ylabel!(plt, latexify("$iv"))
    end

    if !haskey(kw, :title)
        time = round(snap["time"], digits=2)
        title!(plt, "time = $time | pt = $pt")
    end
    
    return plt
end


"""
    histogram_along(snap::Dict, pt::Array=[0.5, 0.5, 0.5]; iv = 0, dir = 1, verbose = 0, kw...)

Plot a histogram of quantity `iv` along direction `dir` through point `pt`. Accepts any keyword argument
accepted by the `Plots.histogram` function.

# Examples
```
julia> histogram_along(snap, [13.4, 20.5, -15.0], dir=3, iv="d", norm=:pdf, label="density") 
```
"""
function histogram_along(snap::Dict, pt=[0.5, 0.5, 0.5]; iv = 0, dir = 1, verbose = 0, kw...)

    kw = Dict(kw)

    kv = Dict{Symbol, Any}(:verbose => 0,
                           :all => false, :iv => 0,
                           :i4 => 0, :label => nothing)

    _kw_extract(kw, kv)
    pt = Float64.(pt)
    data = values_along(snap, pt, iv=iv, dir=dir, all=kv[:all], verbose=verbose)


    if isnothing(kv[:label])
        kv[:label] = latexify("$iv")
    end
    hist = histogram(data, label=kv[:label]; kw...)


    if !haskey(kw, :ylabel)
        if !haskey(kw, :norm)
            ylabel!(hist, "Counts")
        elseif kw[:norm] == :pdf
            ylabel!(hist, "probability density")
        elseif kw[:norm] == :probability
            ylabel!(hist, "probability")
        end
    end

    if !haskey(kw, :title)
        time = round(snap["time"], digits=2)
        title!(hist, "time = $time | pt = $pt")
    end

    return hist
end

"""
    sliceplot(snap::Dict; x = nothing, y = nothing, z = nothing, unigrid = true, kw...)


Make a sliceplot of quantity `iv` at the given coordinate. Accepts any keyword argument
supported by the `Plots.plot()` function, in addition to:

# Arguments
- `snap::Dict`: snapshot
- `iv::Union{Int, String, Char}`: what quantity to plot. Default `0`.
- `grids::Bool`: whether to show patch grids. Default `false`.
- `width::Tuple{Int, Int}`: width of the sub-domain axes. Default `nothing` (whole plane is plotted). 
- `center::Tuple{Int, Int}`:  where to center the plot. Default `nothing` (center is not moved)
- `dims::Union{Tuple{Int, Int}, Int}`: size of interpolated slice if experiment uses adaptive mesh-refinement. 
                                       Can be either tuple indicating size in each dimension, or
                                       an integer indicating same size for both dimensions. Default `nothing`.


# Examples
```
julia> sliceplot(snap, iv="ekin", z=-15.0) # simple sliceplot at z=-15 of kinetic energy
julia> sliceplot(snap, iv="sqrt(bx^2 + by^2 + bz^2)", x=10) # also accepts expressions 
julia> sliceplot(snap, iv="d", z=1.0, center=(5, 5), width=(4, 2)) # re-centre and zoom in
julia> sliceplot(snap, iv="d", z=1.0, resample=true, dims=(400, 600)) # resize data to size (400, 600)
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
                           :width => nothing, :dims => nothing,
                           :center => nothing, :log =>  x -> x,
                           :span => nothing)


    _kw_extract(kw, kv)

    
    iv = kv[:iv]
    verbose = kv[:verbose]

    # get the axis where to slice and the resulting plane
    xyz = [x, y, z]
    dirs = [(x, "x", 1), (y, "y", 2), (z, "z", 3)]
    axis = getindex(dirs, xyz .!= nothing)[1]   # axis at which we are slicing
    planeDirs = getindex(dirs, xyz .== nothing) # axes of plane
    ax1idx, ax2idx = planeDirs[1][3], planeDirs[2][3]
    verbose == 1 && @info ("Sliceplot at $(axis[2]) = $(axis[1]) in plane $((planeDirs[1][2], planeDirs[2][2]))")
    
    origin = copy(snap["cartesian"]["origin"])
    Size = copy(snap["cartesian"]["size"])
    deleteat!(origin, axis[3])
    # # deleteat!(Size, axis[3])

    # check if new width is given
    wflag = false
    if !isnothing(kv[:width])
        width = kv[:width]
        verbose == 1 && println("New width: $width")
        wflag = true
    else
        width = [Size[ax1idx], Size[ax2idx]] / 2
    end

    # check if new center is given
    cflag = false
    if !isnothing(kv[:center])
        center = kv[:center]
        verbose == 1 && println("New center: $center")
        cflag = true
    else
        center = origin .+ width
    end

    # realign data if new width or center is given
    if wflag || cflag
        span = ((center[1] - width[1], center[1] + width[1]), (center[2] - width[2], center[2] + width[2]))
    elseif !isnothing(kv[:span]) 
        span = kv[:span]
    else 
        span = nothing
    end

    
    if unigrid && isnothing(kv[:dims])
        data = unigrid_volume(snap, iv=iv, span=span, verbose=verbose)
        verbose > 1 && @info ("Unigrid data with shape $(size(data))")
        d1 = range(center[1]-width[1], center[1]+width[1], length=size(data, 2)) # dimension 1
        d2 = range(center[2]-width[2], center[2]+width[2], length=size(data, 1)) # dimension 2
    else

        d1, d2, data = amr_volume(snap, iv=iv, span=span, dims=kv[:dims], verbose=verbose)
        verbose > 1 && @info ("Mesh refined data with shape $(size(data))")
    end

    data = @. kv[:log](data)

    # set colorbar title if not given as a keyword arguments
    if !haskey(kw, :cbar_title)
        try
            cbar_title = latexify("$iv") # latexify only works with specific math equations
        catch
            cbar_title = L"%$iv" # produces a latexstring, which should work regardless
        end
    else
        cbar_title = L"%$kw[:cbar_title]"
    end

    # set axis ticks (position, label)
    xticks = range(1, stop=length(d1), step=round(Int, length(d1)/8))
    xticks = (xticks |> collect, ceil.(d1[xticks], digits=1)) 
    yticks = range(1, stop=length(d2), step=round(Int, length(d2)/8))
    yticks = (yticks |> collect, ceil.(d2[yticks], digits=1))

    if !haskey(kw, :aspect_ratio)
        try
            kw[:aspect_ratio] = size(data, 2)/size(data, 1)
        catch
            kw[:aspect_ratio] = 1.0
        end
    end
    
    hm = plot(data, xticks=xticks, yticks=yticks, cbar_title=cbar_title; kw...)

    # add grids
    if kv[:grids]
        i = axis[3]
        patches = patches_in(snap, x=x, y=y, z=z)

        for p in patches
            xidx = ax1idx
            yidx = ax2idx

            e = p["extent"][i,:]

            x = abs.([e[1], e[2], e[2], e[2], e[2], e[1], e[1], e[1]] ./ Size[xidx] .* size(data, 2))
            y = abs.([e[3], e[3], e[3], e[4], e[4], e[4], e[4], e[3]] ./ Size[yidx] .* size(data, 1))

            if i !== 2
                plot!(hm, x, y, color=:gray, label=false)
            else
                plot!(hm, y, x, color=:gray, label=false)
            end
        end
    end

    
    # set labels and title if not given as a keyword argument
    if !haskey(kw, :xlabel)
        xlabel = planeDirs[1][2]
        xlabel!(hm, xlabel)
    end

    if !haskey(kw, :ylabel)
        ylabel = planeDirs[2][2]
        ylabel!(hm, ylabel)
    end

    if !haskey(kw, :title)
        pos = (axis[2], axis[1])
        time = round(snap["time"], digits=2)
        try
            title!(latexify("$(pos[1]) = $(pos[2]), t = $time"))
        catch
            title!(L"$(pos[1]) = $(pos[2]), t = $time")
        end
    end


    return hm

end



""" 
    anim_pane(snap; ax=1, nframes=10, unigrid=true, iv=iv, reverse=false
              start=nothing, stop=nothing, verbose=0, savepath = nothing, kw...)
              
Plot and animate a pane through axis `ax` from `start` to `stop` with `nframes` frames of quantity `iv`. Output is saved
to `savepath`. Will go from low to high axis values unless `reverse` is set to true.

#Examples
```
julia> anim_pane(snap, ax=3, iv="d", nframes=100, start=0.0, stop=-10.0, savepath="test.gif", reverse=true)
```
"""
function anim_pane(snap; ax=1, nframes=10, unigrid=true, iv=0, reverse=false,
                   start=nothing, stop=nothing, span=nothing, verbose=0, savepath = nothing, time=5, kw...)

    if isnothing(savepath)
        throw(ArgumentError("savepath can not be nothing. Aborting."))
    end

    if !isnothing(start) && isnothing(stop)
        axvals = LinRange(start, 
                          snap["cartesian"]["origin"][ax]+snap["cartesian"]["size"][ax], 
                          nframes)
    elseif !isnothing(stop) && isnothing(start)
        axvals = LinRange(snap["cartesian"]["origin"][ax], 
                          stop, 
                          nframes)
    elseif !isnothing(stop) && !isnothing(start)
        axvals = LinRange(start, stop, nframes)   
    else
        axvals = LinRange(snap["cartesian"]["origin"][ax], 
                          snap["cartesian"]["origin"][ax]+snap["cartesian"]["size"][ax], 
                          nframes)    
    end

    if reverse
        axvals = Base.reverse(axvals)
    end
    
    verbose == 1 && println("panning through $ax from $(axvals[1]) to $(axvals[end])")

    data = amr_volume(snap, iv=iv, dims=50)
    vmin, vmax = minimum(data), maximum(data)
    data = nothing

    if ax == 1
        anim_pane_x(snap, axvals, iv = iv, nframes = nframes, unigrid = unigrid, savepath=savepath, time=time, span=span, clims=(vmin, vmax); kw...)
    elseif ax == 2
        anim_pane_y(snap, axvals, iv = iv, nframes = nframes, unigrid = unigrid, savepath=savepath, time=time, span=span, clims=(vmin, vmax); kw...)
    else
        anim_pane_z(snap, axvals, iv = iv, nframes = nframes, unigrid = unigrid, savepath=savepath, time=time, span=span, clims=(vmin, vmax); kw...)
    end
end

""" Wrapper for pane animation in x-direction """
function anim_pane_x(snap, x; iv = 0, nframes = 10, verbose = 0, 
                    unigrid=true, savepath=nothing, time=5, clims=nothing, span=nothing, kw...)

    anim = @animate for i = ProgressBar(1:length(x))
        sliceplot(snap, iv=iv, unigrid=unigrid, x=x[i], title = "x = $(round(x[i], digits=2))",
                  verbose=verbose, clims=clims, span=span; kw...)
    end

    gif(anim, savepath, fps=nframes/time)

end

""" Wrapper for pane animation in y-direction """
function anim_pane_y(snap, y; iv = 0, nframes = 10, verbose = 0, 
                     unigrid=true, savepath=nothing, time=5, clims=nothing, span=nothing, kw...)

    anim = @animate for i = ProgressBar(1:length(y))
        sliceplot(snap, iv=iv, unigrid=unigrid, y=y[i], title = "y = $(round(y[i], digits=2))",
                  verbose=verbose, clims=clims, span=span; kw...)
    end

    gif(anim, savepath, fps=nframes/time)
                
end

""" Wrapper for pane animation in z-direction """
function anim_pane_z(snap, z; iv = 0, nframes = 10, verbose = 0,
                     unigrid=true, savepath=nothing, time=5, clims=nothing, span=nothing, kw...)

    anim = @animate for i = ProgressBar(1:length(z))
        sliceplot(snap, iv=iv, unigrid=unigrid, z=z[i], title = "z = $(round(z[i], digits=2))",
                  verbose=verbose, clims=clims, span=span; kw...)
    end

    gif(anim, savepath, fps=nframes/time)
                
end



function volume(snap::Dict; iv = 0, unigrid=true, kw...)
    kw = Dict{Symbol, Any}(kw)
    if !haskey(kw, :linetype) 
        kw[:linetype] = :heatmap 
    end


    kv = Dict{Symbol, Any}(:verbose => 0, :iv => 0,
                           :grids => false,
                           :width => nothing, :dims => nothing,
                           :center => nothing, :log =>  x -> x,
                           :span => nothing)


    _kw_extract(kw, kv)

    
    iv = kv[:iv]
    verbose = kv[:verbose]

    # get the axis where to slice and the resulting plane
    xyz = [x, y, z]
    dirs = [(x, "x", 1), (y, "y", 2), (z, "z", 3)]
    axis = getindex(dirs, xyz .!= nothing)[1]   # axis at which we are slicing
    planeDirs = getindex(dirs, xyz .== nothing) # axes of plane
    ax1idx, ax2idx = planeDirs[1][3], planeDirs[2][3]
    verbose == 1 && @info ("Sliceplot at $(axis[2]) = $(axis[1]) in plane $((planeDirs[1][2], planeDirs[2][2]))")
    
    origin = copy(snap["cartesian"]["origin"])
    Size = copy(snap["cartesian"]["size"])
    deleteat!(origin, axis[3])
    # # deleteat!(Size, axis[3])

    # check if new width is given
    wflag = false
    if !isnothing(kv[:width])
        width = kv[:width]
        verbose == 1 && println("New width: $width")
        wflag = true
    else
        width = [Size[ax1idx], Size[ax2idx]] / 2
    end

    # check if new center is given
    cflag = false
    if !isnothing(kv[:center])
        center = kv[:center]
        verbose == 1 && println("New center: $center")
        cflag = true
    else
        center = origin .+ width
    end

    # realign data if new width or center is given
    if wflag || cflag
        span = ((center[1] - width[1], center[1] + width[1]), (center[2] - width[2], center[2] + width[2]))
    elseif !isnothing(kv[:span]) 
        span = kv[:span]
    else 
        span = nothing
    end

    
    if unigrid && isnothing(kv[:dims])
        data = unigrid_volume(snap, iv=iv, x=x, y=y, z=z, span=span, verbose=verbose)
        verbose > 1 && @info ("Unigrid data with shape $(size(data))")
        d1 = range(center[1]-width[1], center[1]+width[1], length=size(data, 2)) # dimension 1
        d2 = range(center[2]-width[2], center[2]+width[2], length=size(data, 1)) # dimension 2
    else

        d1, d2, data = amr_plane(snap, iv=iv, x=x, y=y, z=z, span=span, dims=kv[:dims], verbose=verbose)
        verbose > 1 && @info ("Mesh refined data with shape $(size(data))")
    end

    data = @. kv[:log](data)

    # set colorbar title if not given as a keyword arguments
    if !haskey(kw, :cbar_title)
        try
            cbar_title = latexify("$iv") # latexify only works with specific math equations
        catch
            cbar_title = L"%$iv" # produces a latexstring, which should work regardless
        end
    else
        cbar_title = L"%$kw[:cbar_title]"
    end

    # set axis ticks (position, label)
    xticks = range(1, stop=length(d1), step=round(Int, length(d1)/8))
    xticks = (xticks |> collect, ceil.(d1[xticks], digits=1)) 
    yticks = range(1, stop=length(d2), step=round(Int, length(d2)/8))
    yticks = (yticks |> collect, ceil.(d2[yticks], digits=1))

    if !haskey(kw, :aspect_ratio)
        try
            kw[:aspect_ratio] = size(data, 2)/size(data, 1)
        catch
            kw[:aspect_ratio] = 1.0
        end
    end
    
    hm = plot(data, xticks=xticks, yticks=yticks, cbar_title=cbar_title; kw...)

    # add grids
    if kv[:grids]
        i = axis[3]
        patches = patches_in(snap, x=x, y=y, z=z)

        for p in patches
            xidx = ax1idx
            yidx = ax2idx

            e = p["extent"][i,:]

            x = abs.([e[1], e[2], e[2], e[2], e[2], e[1], e[1], e[1]] ./ Size[xidx] .* size(data, 2))
            y = abs.([e[3], e[3], e[3], e[4], e[4], e[4], e[4], e[3]] ./ Size[yidx] .* size(data, 1))

            if i !== 2
                plot!(hm, x, y, color=:gray, label=false)
            else
                plot!(hm, y, x, color=:gray, label=false)
            end
        end
    end

    
    # set labels and title if not given as a keyword argument
    if !haskey(kw, :xlabel)
        xlabel = planeDirs[1][2]
        xlabel!(hm, xlabel)
    end

    if !haskey(kw, :ylabel)
        ylabel = planeDirs[2][2]
        ylabel!(hm, ylabel)
    end

    if !haskey(kw, :title)
        pos = (axis[2], axis[1])
        time = round(snap["time"], digits=2)
        try
            title!(latexify("$(pos[1]) = $(pos[2]), t = $time"))
        catch
            title!(L"$(pos[1]) = $(pos[2]), t = $time")
        end
    end


    return hm

"""
    anim_plane(;data="../data", run="", x = nothing, y = nothing, z = nothing, iv=0, 
                tspan=nothing, unigrid=true, step=1, savepath=nothing, verbose = 0, kw...)

Animate a plane of quantity `iv` at a slice `x/y/z`. If `tspan` tuple is not given, all snapshots will be used, otherwise only snapshots
with time within the given `tspan` will be loaded. A frame will be recorded every `step`, and the resulting animation
will be saved to `savepath`. 
"""
function anim_plane(;data="../data", run="", x = nothing, y = nothing, z = nothing, iv=0, 
                    tspan=nothing, unigrid=true, step=1, savepath=nothing, verbose = 0, kw...)

    if isnothing(savepath)
        throw(ArgumentError("savepath can not be nothing. Aborting."))
    end

    snapIDs = sort(get_snapshot_ids(data=data, run=run)) # all snapshot IDs

    # if tspan is given, find the snapshots where time=t0 and time=tmax
    if !isnothing(tspan)
        t0, tmax = tspan
        start, stop = nothing, nothing

        times = get_snapshot_time.(snapIDs, data=data, run=run)
        start = snapIDs[argmin(abs.(times .- t0))]
        stop = snapIDs[argmin(abs.(times .- tmax))]
    else # otherwise use all snapshots
        start = snapIDs[1]
        stop = snapIDs[end]
        t0 = get_snap_time(start, data=data, run=run)
        tmax = get_snap_time(stop, data=data, run=run)
    end
    
    snapIDs = snapIDs[start:step:stop]

    anim = @animate for i in snapIDs
        snap = snapshot(i, data=data, run=run, verbose = verbose)
        sliceplot(snap, iv=iv, unigrid=unigrid, x=x, y=y, z=z,
                  verbose=verbose; kw...)
    end

    gif(anim, savepath)
    
end
