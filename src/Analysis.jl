module Analysis

include("analysis/analysis.jl")
export power_spectrum2d, density2d, horizontal_average, stacked_density, LaplacianMatrix

include("tools/vector_calculus.jl")
export gradient, divergence, curl

end
