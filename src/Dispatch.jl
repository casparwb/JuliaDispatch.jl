module Dispatch
    using JLD
    include("dispatch/_dispatch.jl")

    export snapshot, cache_snapshots_live
end