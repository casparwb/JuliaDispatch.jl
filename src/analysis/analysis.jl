# using KernelDensity
using JuliaDispatch.Dispatch
using JuliaDispatch.Utils
using JuliaDispatch.Select
using StaticArrays, ProgressBars
using Base.Threads

function mean(x)
    return sum(x)/length(x)
end

"""
    average(snap::Dict; iv::Union{Int, String})

Average of quantity `iv` over all patches
"""
function average(snap; iv=0, all=false)

    data = 0.0
    vol = 0

    n = all ? "gn" : "n"
    if nthreads() > 1
        mydata = zeros(Float64, nthreads())
        myvol = zeros(Int, nthreads())
        @threads for patch in snap["patches"]
            mydata[threadid()] += sum(patch["var"](iv, all=all))# |> sum
            myvol[threadid()] += patch[n] |> prod
        end
        return sum(mydata)/sum(myvol)

    else
        data = 0.0
        vol = 0
        for patch in snap["patches"]
            data += sum(patch["var"](iv, all=all))# |> sum
            vol += patch[n] |> prod
        end
        return data/vol
    end
end


"""
    time_evolution_average(;iv = 0, run="", data="../data",
                            all=false, verbose=0, tspan=nothing)

Time evolution of quantity `iv` averaged over entire computational domain.
"""
function time_evolution_average(;iv = 0, run="", data="../data",
                                 all=false, verbose=0, tspan=nothing)

    verbose == 1 && @info "Time evolution average of quantity $iv"
    times = get_snapshot_time(data=data, run=run, tspan=tspan)
    snapIDs = get_snapshot_ids(data=data, run=run, tspan=tspan)

    averages = Vector{Float64}(undef, length(times))
    for index ∈ ProgressBar(eachindex(times))
        averages[index] = average(snapshot(snapIDs[index], data=data, run=run, progress=false, suppress=true), iv=iv, all=all)
    end

    return times, averages
end

"""
    time_evolution_plaverage(;iv = 0, run="", data="../data", dir=3, tspan=nothing,
               all=false, verbose=0, nslices=nothing)

Time evolution of plane average of quantity `iv` along direction `dir`. Returns 
a `Vector` containing the times of each snapshot, a `Vector` containing the axis coordinates along direction
`dir`, and a `Matrix` with the averages at each coordinate (axis 1) at each time (axis 2). 
"""
function time_evolution_plaverage(;iv = 0, run="", data="../data", dir=3, tspan=nothing,
                                   all=false, verbose=0, nslices=nothing)

    times = get_snapshot_time(data=data, run=run, tspan=tspan)
    snapIDs = get_snapshot_ids(data=data, run=run, tspan=tspan)

    verbose == 1 && @info "Time evolution plane average in direction $dir from t=$(times[1]) to t=$(times[end])"

    averages = Vector{Vector{Float32}}(undef, length(times))
    axes = nothing
    for index ∈ ProgressBar(eachindex(times))
        snap = snapshot(snapIDs[index], data=data, run=run, progress=false, verbose=-1)
        axes, averages[index] = plaverage(snap, iv=iv, dir=dir, all=all, verbose=verbose, nslices=nslices)
    end
    
    try
        averages = hcat(averages...)
    catch e
        @warn "Could not concatenate array. Returning array of arrays"
    end

    return axes, times, averages
end

"""
    plaverage(snap; iv=0, dir=3, all=false, nslices=nothing)

Compute and return the average value of quantity `iv` in each plane slice perpendicular
to direction `dir`.
"""
function plaverage(snap; iv=0, dir=3, all=false, verbose=0, nslices=nothing)

    if isnothing(nslices)
        try
            axis_length = snap["datashape"][dir]
        catch e
            @error "Snap has no datashape property. Please give nslices keyword argument value."
            return nothing
        end
    else
        axis_length = nslices
    end

    verbose == 1 && @info "Plane average in direction $dir with $axis_length slices"

    axis = range(snap["cartesian"]["origin"][dir], 
                 snap["cartesian"]["origin"][dir] + snap["cartesian"]["size"][dir],
                 length=axis_length)

    quantity_average = zeros(Float64, length(axis))
    
    for index ∈ eachindex(axis)
        x, y, z = (i == dir ? (axis[index]) : nothing for i = 1:3)
        patches = patches_in(snap, x=x, y=y, z=z)
        Base.Threads.@threads for patch in patches
            data = plane(patch, iv=iv, x=x, y=y, z=z, all=all)
            quantity_average[index] += sum(data)/length(data)
        end
    end

    return axis, quantity_average
end

"""
    stacked_density(snap::Dict, N::Int, vmin, vmax; iv=0, dir = 3, all=false, nbins=50)

Return a stacked kernel density estimate of quantity `iv` in each plane perpendicular to axis `dir`.

WORK IN PROGRESS! 
"""
function stacked_density(snap::Dict, N::Int, vmin, vmax; iv=0, dir = 3, all=false, nbins=50)

    Size = snap["cartesian"]["size"]
    origin = snap["cartesian"]["origin"]

    ydata = range(origin[dir], origin[dir]+Size[dir], length=N+1)
    ydata = ydata[1:end-1]

    if nbins < 100
        xdata = SVector{nbins}(collect(range(0.8*vmin, 1.2*vmax, length=nbins)))
        bins = SVector{nbins-1}([(xdata[i], xdata[i+1]) for i = 1:nbins-1])
    else
        xdata = range(0.8*vmin, 1.2*vmax, length=nbins)
        bins = [(xdata[i], xdata[i+1]) for i = 1:nbins-1]
    end

    density = zeros(Float64, N, nbins)

    @inbounds for idx = 1:N
        x, y, z = [i == dir ? ydata[idx] : nothing for i = 1:3]
        patches = patches_in(snap, x = x, y = y, z = z)
        @threads for patch in patches
            slice = plane(patch, iv=iv, z=ydata[idx])
            @inbounds for col in eachcol(slice)
                KDE = kde(col)
                var = KDE.x
                dens = KDE.density

                @inbounds for (i, bin) in enumerate(bins)
                    idxs = findall(bin[1] .<= var .<= bin[2])
                    isempty(idxs) && continue
                    density[idx, i] += sum(dens[idxs])
                end
            end
        end
    end

    return xdata, ydata, density
end
