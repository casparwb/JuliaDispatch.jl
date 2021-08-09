# using FFTW, KernelDensity
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

function time_evolution_plaverage(;iv = 0, run="", data="../data", dir=3, tspan=nothing,
                                   all=false, verbose=0, nslices=nothing)

    times = get_snapshot_time(data=data, run=run, tspan=tspan)
    snapIDs = get_snapshot_ids(data=data, run=run, tspan=tspan)

    verbose == 1 && @info "Time evolution plane average in direction $dir from t=$(times[1]) to t=$(times[end])"

    averages = Vector{Vector{Float32}}(undef, length(times))
    axes = similar(averages)
    for index ∈ ProgressBar(eachindex(times))
        snap = snapshot(snapIDs[index], data=data, run=run, progress=false, suppress=true)
        axes[index], averages[index] = plaverage(snap, iv=iv, dir=dir, all=all, verbose=verbose, nslices=nslices)
    end

    if Base.all(length.(axes) .== length(axes[1])) 
        axes = axes[1] 
        averages_matrix = Matrix{Float32}(undef, length(axes), length(times))
        for index in eachindex(times)
            averages_matrix[:,index] = averages[index]
        end

        return axes, times, averages_matrix
    end

    return times, axes, averages
end

"""
    plaverage(snap; iv=0, dir=3, all=false, nslices=nothing)

Compute and return the average value of quantity `iv` in each plane slice perpendicular
to direction `dir`.

"""
function plaverage(snap; iv=0, dir=3, all=false, verbose=0, nslices=nothing)

    if isnothing(nslices)
        axis_length = snap["datashape"][dir]
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
    stacked_density(snap::Dict, N::Int; iv::Union{Int, String}, dir::Int,
                    all::Bool, nbins::Int, Log::Bool)

Return a stacked kernel density estimate of quantity `iv` in each plane perpendicular to axis `dir`.

# Arguments:
- `snap::Dict`: snapshot object
- `N::Int`: number of plane slices

# Kwargs:
- `iv::Union{Int, String}`: variable name, default 0
- `dir::Int`: axis at which to slice data, default 3 (z)
- `nbins::Int`: number of histogram bins, default 50
- `Log::Bool`: whether to return log of x-axis, default false
- `all::Bool`: whether to include guard zones, default false

# Returns:
- Array{Float, 1} of length nbins, binned data values in ascending order
- Array{Float, 1} of length N, axis values for axis dir
- Array{Float, 2} of size (N, nbins), kernel density at each plane
"""
function stacked_density(snap::Dict, N::Int; iv=0, dir = 3, all=false, nbins=50,
                         Log = false)

    Size = snap["cartesian"]["size"]
    origin = snap["cartesian"]["origin"]

    ydata = range(origin[dir], origin[dir]+Size[dir], length=N+1)
    ydata = ydata[1:end-1]
    minVal, maxVal = Inf, -Inf

    @inbounds for patch in snap["patches"]
        data = box(patch, iv=iv)
        Min, Max = minimum(data), maximum(data)

        if Min < minVal
            minVal = Min
        end

        if Max > maxVal
            maxVal = Max
        end
    end

    if nbins < 100
        xdata = SVector{nbins}(collect(range(0.8*minVal, 1.2*maxVal, length=nbins)))
        bins = SVector{nbins-1}([(xdata[i], xdata[i+1]) for i = 1:nbins-1])
    else
        xdata = range(0.8*minVal, 1.2*maxVal, length=nbins)
        bins = [(xdata[i], xdata[i+1]) for i = 1:nbins-1]
    end

    density = zeros(Float64, N, nbins)

    @inbounds for idx = 1:N
        patches = patches_in(snap, z=ydata[idx])
        @inbounds for patch in patches
            slice = plane(patch, iv=iv, z=ydata[idx])
            @inbounds for col in eachcol(slice)
                KDE = kde(col)
                var = KDE.x
                dens = KDE.density

                for (i, bin) in enumerate(bins)
                    idxs = findall(bin[1] .<= var .<= bin[2])
                    isempty(idxs) && continue
                    density[idx, i] += sum(dens[idxs])
                end
            end
        end
    end

    Log ? xdata = log10.(xdata) : nothing
    return xdata, ydata, log10.(density)
end


"""
    stacked_density(data::Array{Float, 3}; dir::Int, nbins::int, Log::Bool, step::Int)

Return a stacked kernel density estimate of each plane perpendicular to
axis `dir`.

Arguments:
-----------
    - data: Array{Float, 3}, 3d array with data

Kwargs:
-----------
    - dir:   Int, axis at which to slice data, default 3 (z)
    - nbins: Int, number of histogram bins, default 50
    - Log:   Bool, whether to return log of x-axis, default false
    - step:  Int, step size for slices along dir, default 1

Returns:
----------
    - Array{Float, 1} of length nbins, binned data values in ascending order
    - Array{Float, 2} of size (size(data[dir]), nbins), kernel density at each plane
"""
function stacked_density(data::Array{T, 3} where T <: Number;
                         dir = 3, nbins=50, Log=true, step=1)



    sze = size(data)[dir]
    xdata = range(minimum(data)*0.8, maximum(data)*1.2, length=nbins)
    bins = [(xdata[i], xdata[i+1]) for i = 1:nbins-1]
    density = zeros(sze, nbins)

    for idx = 1:step:sze
        slice = view(data, :,:,idx)
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

    if Log
        return log10.(xdata), log10.(density)
    else
        return xdata, log10.(density)
    end
end
