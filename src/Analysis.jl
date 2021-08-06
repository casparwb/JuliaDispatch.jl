module Analysis

include("analysis/analysis.jl")
export power_spectrum2d, density2d, horizontal_average, stacked_density, average, 
       plaverage, time_evolution_plaverage, time_evolution_average

include("tools/vector_calculus.jl")
export gradient, divergence, curl

end
