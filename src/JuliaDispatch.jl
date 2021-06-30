module JuliaDispatch

    include("_dispatch.jl")
    include("graphics/_graphics.jl")
    include("select/buffers.jl")
    include("analysis.jl")
    include("tools/vector_calculus.jl")

    export unigrid_volume, unigrid_plane, amr_volume, amr_plane, resample
    export snapshot, sliceplot, plot_values_along
    export power_spectrum2d, density2d, horizontal_average, stacked_density
    export gradient, divergence, curl
end # end module