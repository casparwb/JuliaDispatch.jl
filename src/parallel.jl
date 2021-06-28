function stacked_density_parallel(axvals::AbstractArray{T, 1} where T,
                                  data::Array{T, 3} where T;
                                  dir = 3, verbose=0, nbins=10, _log=true)



    xdata = range(minimum(data)*0.5, maximum(data)*1.5, length=nbins)
    ydata = axvals
    bins = [(xdata[i], xdata[i+1]) for i = 1:nbins-1]
    densData = zeros(length(ydata), nbins)

    shp = size(data)
    dirSize = shp[dir]

    @assert length(ydata) == dirSize "axvals must be same length as size(data) in dimension $dir"

    @distributed for idx = 1:dirSize
        Plane = data[:,:,idx]'

        @inbounds for j = 1:size(Plane)[1]
            # var, dens = xy(ash(Plane[j,:]))
            KDE = kde(Plane[j,:])
            var = KDE.x
            dens = KDE.density
            # var = Plane[j,:]

            @inbounds for (i, bin) in enumerate(bins)
                idxs = findall(bin[1] .<= var .<= bin[2])
                isempty(idxs) && continue
                # densData[idx, i] += length(idxs)
                densData[idx, i] += sum(dens[idxs])
            end
        end
    end

    if _log
        return log10.(xdata), ydata, log10.(densData)
    else
        return xdata, ydata, log10.(denData)
    end
end
