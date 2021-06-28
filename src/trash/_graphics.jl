# using Plots, PlotThemes
# theme(:juno)
using Makie


include("../select/_select.jl")

function plot_values_along(pp, pt=[0.5, 0.5, 0.5]; kw...)
    """ Plot values along direction dir={1,2,3}, through point pt=[x,y,z] """
    kw = Dict(kw)
    kv = Dict{Symbol, Union{String, Int, Bool}}(:dir => 1, :verbose => 0,
                                                :all => false, :iv => 0,
                                                :i4 => 0, :var => nothing)


    _kw_extract(kw, kv)

    data = values_along(pp, pt, iv=kv[:iv], dir=kv[:dir], all=kv[:all])

    plt = plot(data, label="$iv")
    if kv[:dir] == 1
        dir = "x"
        xlabel!("x")
        ylabel!("z")
    elseif kv[:dir] == 2
        dir = "y"
        xlabel!("y")
        ylabel!("z")
    elseif kv[:dir] == 3
        dir = "z"
        xlabel!("x")
        ylabel!("y")
    end

    title!("$iv values along $dir")

    plt

end

function plot_values_along!(pp, pt=[0.5, 0.5, 0.5]; kw...)
    """ Plot values along direction dir={0,1,2}, through point pt=[x,y,z] """
    kw = Dict(kw)
    kv = Dict(:dir => 1, :verbose => 0, :all => false,
              :iv => 0, :i4 => 0, :var => nothing)


    _kw_extract(kw, kv)

    data = values_along(pp, pt, iv=kv[:iv], dir=kv[:dir], all=kv[:all])

    Plots.plot!(data, label=false)

end

function _kw_extract(kw, dict)
    """ if key from dict occur in kw, pop them """
        for (key, value) in dict
            dict[key] = pop!(kw, key, value)
            try
                println(kw[key])
            catch
                nothing
            end
        end
    # return kw, dict
end

# function plot(f; kw...)
#     """ plot(f) allows f to be (x, y) tuple """
#
#     if typeof(f) <: Tuple
#         x, y = f
#         """ PLOT """ # plt.plot(x, y, kw...)
#     else
#         # plt.plot(f,kw...)
#     end
# end


function get_n_snapshots(run="", data="../data")
    """
    Returns the total number of snapshots in the given data folder
    """

    datadir = joinpath(run, data)
    nsnaps = length([dir for dir in readdir(datadir) if startswith(dir, "0")])
end


function image_plane(snap; x = nothing, y = nothing, z = 0.5, kw...)

     kw = Dict(kw)
     kv = Dict{Symbol, Union{String, Int, Bool}}(:verbose => 0, :iv => 0,
                                                 :grids => false, :cbar => false,
                                                 :title => 0)


     _kw_extract(kw, kv)

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

    fmin < 0 ? cmap = :coolwarm : cmap = :viridis

    sortedkeys = sort(collect(keys(ll)))

    data = reshape(zeros(Float32, datasize), Int(sqrt(datasize)), :)
    datashp = size(data)

    for id in sortedkeys
        im, e = ll[id][1:2]

        xs = Int.(e[1:2] .* datashp)
        ys = Int.(e[3:4] .* datashp)

        data[(xs[1] + 1):xs[end], (ys[1] + 1):ys[end]] = im
    end

    xx = collect(LinRange(xmin, xmax, datashp[1]))
    yy = collect(LinRange(ymin, ymax, datashp[2]))



    scene = heatmap(xx, yy, data, interpolate=interpolate, colormap=cmap)
    # xlims!(scene, emin[1], emax[1])
    # ylims!(scene, emin[2], emax[2])

    if title_ != nothing
        title(scene, title)
    else
        title(scene, "$(sizedir[i]) = $(xyz[i])")
    end

    cs = colorlegend(scene[end], range=[fmin, fmax])#, label="$iv")

    if grids
        for p in patches
            e = p["extent"][i,:]

            x = [e[1], e[1], e[2], e[2], e[1]]
            y  =[e[3], e[4], e[4], e[3], e[3]]

            plot!(scene, x, y, color=:white)
        end
    end

    xlabel!(scene, labels[i][1])
    ylabel!(scene, labels[i][2])

    return vbox(scene, cs)
    

end


# function amr_plane(snap; iv="d", z=0.5, log=false, mesh=false, ident=false,
#                    cmap=nothing, title=nothing, zero=false, vmin=nothing,
#                    vmax=nothing, verbose = 0)
#
#     scene = Makie.scene()
#     contains(p, z) = p["level"] >= 6 || abs(z - p["positions"][3]) <= p["size"][3]/2
#
#     function outline!(scene, e)
#         """ Draw outline of patch with extent e """
#         x = [e[1],e[1],e[2],e[2],e[1]]
#         y = [e[3],e[4],e[4],e[3],e[3]]
#         plot!(scene, x,y,color=:gray)
#     end
#
#     function interpolate(patch, iv, z, log=false)
#         n = p["n"][3]
#         dz = p["size"][3]/n
#         z0 = p["position"][3] - p["size"][3]/2 + dz/2
#         w = (z - z0)/dz
#
#         i = Int(round(w))
#         i = max(, min(i, n-2))
#         w -= i
#
#         f = p["var"](iv)
#
#         if iv == 0 || iv == 4 || iv == "d" || iv = "e"
#             log = true
#         end
#
#         log ? f = log.(f) : nothing
#
#         i == 0 ? f = (1 - w)*f[:,:,end] + w*f[:,:,i+1] : nothing
#
#         log ? f = log.(f) : nothing
#
#         return f
#     end
#
#     patches = [p for p in snap["patches"] if contains(p, z)]
#
#     vmin1 = 1e20
#     vmax1 = 0.0
#
#     data = []
#     for patch in patches
#         d = interpolate(patch, iv, z)
#         e = p["extent"][3,:]
#         push!(data, (d, e))
#
#         vmin1 = min(vmin1, minimum(d))
#         vmax1 = max(vmin1, maximum(d))
#
#         verbose > 0 ? println("$(p["id"]), $(minimum(d)), $(maximum(d))") : nothing
#     end
#
#     vmin == nothing ? vmin = vmin1 : nothing
#     vmax == nothing ? vmax = vmax1 : nothing
#
#     if zero
#         vmax = max(abs(vmax), abs(vmin))
#         vmin = -vmax
#     end
#
#     xlims!(scene, 0, 1)
#     ylims!(scene, 0, 1)
#
#     if log
#         vmin = log10(vmin)
#         vmax = log10(vmax)
#     end
#
#     for (d, e) in data
#         if log
#             d = log10.(d)
#         end
#
#         heatmap(scene)
#     end
#
# end
