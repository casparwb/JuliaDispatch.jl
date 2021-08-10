

function curl(V, X)

    ∂V₁∂X₂ = gradient(V[1], spacing=X[2], dim=2)
    ∂V₁∂X₃ = gradient(V[1], spacing=X[3], dim=3)

    ∂V₂∂X₃ = gradient(V[2], spacing=X[3], dim=3)
    ∂V₂∂X₁ = gradient(V[2], spacing=X[1], dim=1)

    ∂V₃∂X₂ = gradient(V[3], spacing=X[2], dim=2)
    ∂V₃∂X₁ = gradient(V[3], spacing=X[1], dim=1)

    return [∂V₃∂X₂ .- ∂V₂∂X₃,
            ∂V₁∂X₃ .- ∂V₃∂X₁,
            ∂V₂∂X₁ .- ∂V₁∂X₂]

end


function gradient(F; spacing=1, dim=1)

    not_dims = [size(F)...]
    not_dims[dim] = 1
    grad = cat(zeros(not_dims...), diff(F, dims=dim), dims=dim)

    isa(spacing, Number) && return grad / spacing

    @assert length(spacing) == size(grad, dim) "Given spacing array does not match F along dimension $dim. Got $(length(spacing)), must be $(size(F, dim))"
    
    try
        for index ∈ axes(grad, dim)
            data = selectdim(grad, dim, index)
            data ./= spacing[index]
        end
        return grad
    catch e
        @error "Could not compute gradient. $e"
        return nothing
    end
end

const ∇ = gradient