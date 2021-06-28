@everywhere using KernelDensity
@everywhere using SharedArrays
# @everywhere include("utilities/julia/_dispatch.jl")

function stacked_density_p(axvals::AbstractArray{T, 1} where T,
                           data::Array{T, 3} where T;
                           dir = 3, verbose=0, nbins=10, _log=true)



    xdata = range(minimum(data)*0.5, maximum(data)*1.5, length=nbins)
    ydata = axvals
    bins = [(xdata[i], xdata[i+1]) for i = 1:nbins-1]

    # @everywhere densData = SharedArray{Float64}((300, 20), init=zeros())
    densData = SharedArray{Float64}((length(axvals), nbins), init=zeros())

    @assert length(ydata) == size(data)[dir] "axvals must be same length as size(data) in dimension $dir"

    @sync @distributed for idx = 1:size(data)[dir]
        slice = view(data, :,:,idx)
    @inbounds for col in eachcol(slice)
            KDE = kde(col)
            var = KDE.x
            dens = KDE.density

            @inbounds for (i, bin) in enumerate(bins)
                idxs = findall(bin[1] .<= var .<= bin[2])
                isempty(idxs) && continue
                densData[idx, i] += sum(dens[idxs])
            end
        end
    end

    if _log
        return log10.(xdata), ydata, log10.(densData)
    else
        return xdata, ydata, log10.(densData)
    end
end

@everywhere function box(patch; iv = 0, all=false, verbose = 0)


    if patch["guard_zones"] && !all
        li = patch["li"]  # lower inner
        ui = patch["ui"]  # upper inner
    elseif all || !patch["guard_zones"]
        li = ones(Int, 3)
        ui = patch["gn"]
    end

    f = patch["var"](iv)[li[1]:ui[1], li[2]:ui[2], li[3]:ui[3]]

    return f
end


function stacked_density_p(snap::Dict, N::Int; iv=0, dir = 3, all=false, nbins=50,
                  Log = false)

    Size = snap["cartesian"]["size"]
    origin = snap["cartesian"]["origin"]

    ydata = range(origin[dir], origin[dir]+Size[dir], length=N+1)
    ydata = ydata[1:end-1]

    minmax = SharedArray{Float32}([Inf, -Inf])
    minVal, maxVal = Inf, -Inf
    minVal = @sync @distributed min for patch in snap["patches"]
        data = box(patch, iv=iv)
        # Min, Max = minimum(data), maximum(data)
        # minVal = min(minVal, Min)
        # maxVal = max(maxVal, Max)
        minimum(data)
    end

    println(minVal)
    println(maxVal)

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
