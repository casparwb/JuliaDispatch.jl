module analysis

include("_analysis.jl")
export power_spectrum2d, density2d, horizontal_average, stacked_density

include("vector_calculus.jl")
export gradient, divergence, curl

end
