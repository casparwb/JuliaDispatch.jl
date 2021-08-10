module Analysis

include("analysis/analysis.jl")
export average, plaverage, time_evolution_plaverage, time_evolution_average

include("tools/vector_calculus.jl")
export gradient, diff_quotient, ∇

end
