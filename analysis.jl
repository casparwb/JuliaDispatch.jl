using FFTW, KernelDensity
include("_dispatch.jl")

function power_spectrum2d(data; kw...)
    #kw = Dict(kw)
    ft2 = abs.(fft(data).^2)
    shp = size(data)

    k1 = collect(1:shp[2])' .* ones(shp[1])
    k2 = collect(1:shp[1]) .* ones(shp[2])'

    k = sqrt.(k1.^2 .+ k2.^2)

    a = 2
    k0 = 1/sqrt(a)
    k1 = 1/sqrt(a)

    power = Array{Float32, 1}([])
    kk = similar(power)

    sizehint!(power, shp[1])
    sizehint!(kk, shp[1])

    while k1 <= shp[1]÷2
        push!(kk, sqrt(k0*k1))

        w = findall(k0 .< k .<= k1)
        push!(power, sum(ft2[w]))

        k0 = k1
        k1 *= a
    end

    # plot(kk, power, yscale=:log10, xscale=:log10)
    return kk, power
end

"""
    average(snap::Dict; iv::Union{Int, String})

Average of quantity `iv` over all patches
"""
function average(snap; iv=0)

    f = 0.0
    vol = 0.0
    for patch in snap["patches"]
        data = patch["var"](iv)    # patch data
        dvol = prod(patch["ds"])   # volume of each cell
        vol += dvol*prod(size(data))
        f += dvol*sum(patch["var"](iv, all=false))
    end
    return f/vol
end

function time_average(;iv = 0, run="", data="../data/")

    snapshots = nothing

    t_aver = 0.0
    time = Array{Float32, 1}(undef, nsnaps)

    nsnaps = get_n_snapshots(run=run, data=data)
    for iout = 0:nsnaps-1
        snap = snapshot(iout, run=run, data=data)
        time[iout+1] = snap["time"]
        t_aver += average(snap, iv=iv)
    end

    return t_aver ./ (time[end] - time[1])

end

function time_average(snaps, tlims=nothing; iv=0, run="", data="../data/")

    if tlims != nothing
        snaps = [snap for snap in snaps if (tlims[1] <= snap["time"] <= tlims[2])]
    end

    times = [snap["time"] for snap in snaps]

    tmin = minimum(times)
    tmax = maximum(times)
    dt = tmax - tmin

    t_aver = 0.0
    idx = 1
    while time <= tlims[2] || idx <= lenghth(snaps)
         t_aver += average(snaps[idx])
         idx += 1
    end

     return t_aver ./ dt
end

""" Horizontal average of iv in direction dir over all times """
function time_haverage(;iv = 0, run="", data="../data", dir=3,
                        i4=0, all=false, verbose=0)

    nsnaps = get_n_snapshots(run=run, data=data)

    # axVals = Array{Float32, 1}([])
    times = Array{Float32, 1}(undef, nsnaps)

    snap0 = snapshot(0, run=run, data=data)
    axVals, hav = horizontal_average(snap0, iv=iv, dir=dir, all=all, i4=i4)
    havs = Array{Float32, 2}(undef, nsnaps, length(hav))
    havs[1,:] = hav
    times[1] = snap0["time"]

    for iout ∈ 1:nsnaps-1
        snap = snapshot(iout, run=run, data=data)
        times[iout+1] = snap["time"]
        havs[iout+1,:] = horizontal_average(snap, iv=iv, dir=dir, all=all, i4=i4)[2]
        # push!(havs, hav)
    end

    return times, axVals, havs
end

"""
    haver(snap::Dict; iv::Union{Int, String}, dir::Int, i4::Int, all::Bool, Verbose::Int)

Compute and return the average value of quantity `iv` in each plane slice perpendicular
to direction `dir`.

Arguments:
--------------
    -
"""
function haver(snap; iv=0, dir=3, i4=0, all=false, verbose=0)
    patches = snap["patches"]
    patch0 = patches[1]

    all ? posVals = "xyz" : posVals = "xyzi"

    xx = Array{Float32, 1}([])
    h_av = Array{Float32, 1}([])
    sizehint!(xx, patch0["gn"][dir]*length(patches))
    sizehint!(h_av, patch0["gn"][dir]*length(patches))

    jv = map_var(patch0, iv)
    jv == 0 ? jv += 1 : nothing
    for patch in patches
        # patch["no_mans_land"] ? es = 0.5 : es = 0.0
        # if jv == -1
        #     es += patch["idx"]["h"][dir, end]
        # else
        #     es += patch["idx"]["h"][dir, jv]
        # end
        all ? n = patch["gn"] : n = patch["n"]

        rr = patch[posVals][dir]

        f = nothing
        if iv in patch["all_keys"]
            data = box(patch, iv=iv, all=all)
            for (i, f) in enumerate(eachslice(data, dims=dir))
                push!(xx, rr[i])
                push!(h_av, sum(f)/length(f))
            end
        end
    end

    # sort the results
    sortidxs = sortperm(xx)
    xx = xx[sortidxs]
    h_av = h_av[sortidxs]

    # harvest (remove repeated values)
    ss = Array{Float32, 1}([])
    hh = similar(ss)

    i = 1
    while i <= length(xx)
        x0 = xx[i]
        n = 0
        ha = 0.0
        while i <= length(xx) && xx[i] == x0
            ha += h_av[i]
            i += 1
            n += 1
        end
        push!(ss, x0)
        push!(hh, ha/n)
    end

    return ss, hh
end

"""
    haver(data::Array{Float, 3}; dir::Int)

Return horizontal average of `data` in all planes along direction `dir`.
"""
function haver(data::Array{T, 3} where T <: Number; dir=3)
    aver = [sum(slice)/length(slice) for slice in eachslice(data, dims=dir)]
end

"""
    stacked_density(snap::Dict, N::Int; iv::Union{Int, String}, dir::Int,
                    all::Bool, nbins::Int, Log::Bool)

Return a stacked kernel density estimate of quantity `iv` in each plane perpendicular to axis `dir`.

Arguments:
-----------
    - snap: Dict, snapshot object
    - N:    Int, number of plane slices

Kwargs:
-----------
    - iv:    Int/String, variable name, default 0
    - dir:   Int, axis at which to slice data, default 3 (z)
    - nbins: Int, number of histogram bins, default 50
    - Log:   Bool, whether to return log of x-axis, default false
    - all:   Bool, whether to include guard zones, default false
Returns:
----------
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


mean(A) = sum(A)/length(A)

function density2d(snap; iv=0, ax=1, x=nothing, y=nothing, z=nothing, unigrid=false)
    """
    Compute the mean density of iv along dimension ax at a slice x/y/z
    """
    if unigrid
        Plane = unigrid_plane(snap, iv=iv, x=x, y=y, z=z)
    else
        Plane = amr_plane(snap, iv=iv, x=x, y=y, z=z)'
    end

    dim = nothing
    if ax == 1
        if z != nothing
            dim = 1
        elseif y != nothing
            dim = 1
        end
    elseif ax == 2
        if x != nothing
            dim = 1
        elseif z != nothing
            dim = 2
        end
    elseif ax == 3
        if y != nothing
            dim = 2
        elseif x != nothing
            dim = 2
        end
    end

    start = snap["cartesian"]["origin"][ax]
    xdata = range(start, snap["cartesian"]["size"][ax]+start, length=size(Plane)[dim])

    # ydata = nothing
    # if dim == 1
    #     ydata = [sum(Plane[i,:])/length(Plane[i,:]) for i = 1:size(Plane)[1]]
    # else
    #     ydata = [sum(Plane[:,i])/length(Plane[:,i]) for i = 1:size(Plane)[2]]
    # end
    ydata = [sum(Plane[i,:])/length(Plane[i,:]) for i = 1:size(Plane)[dim]]

    KDE = kde((xdata, ydata))

    return KDE.x, KDE.y, KDE.density
end
