

"""
    gradient(F::AbstractArray; spacing::Union{AbstractArray, Number}=1)

Compute the difference quotient of `F` along all dimensions. `spacing` can be a number
specifying uniform sample distance for all dimensions, or an array or, a collection of 
`ndims(F)` arrays for each value in each dimension.
"""
function gradient(F; spacing=1)

    ndims(F) == 1 && return diff_quotient(F, spacing=spacing, dim=1)
    dims = 1:ndims(F)
    if length(spacing) == ndims(F)
        return [diff_quotient(F, spacing=spacing_, dim=i) for (spacing_, i) in zip(spacing, dims)]
    else
        return [diff_quotient(F, spacing=spacing, dim=i) for i in dims]
    end
end

"""
    diff_quotient(F::AbstractArray; spacing::Union{AbstractArray, Number}=1, dim=1)

Compute the difference quotient of `F` along dimension `dim` using finite difference with `spacing`
giving the coordinate distance between elements in `F`. `spacing` can be a number, which
is equivalent to equal spacing among all points, or it can be an array with spacing between the values.
Must then be equal to `size(F, dim)`.  
"""
function diff_quotient(F; spacing=1, dim=1)
    
    @assert dim <= ndims(F) "dim = $dim was given, but F has only $(ndims(F)) dimensions."
    all_dims = [size(F)...]
    all_dims[dim] = 1
    quotient = cat(zeros(all_dims...), diff(F, dims=dim), dims=dim)

    isa(spacing, Number) && return quotient / spacing

    @assert length(spacing) == size(quotient, dim) "Given spacing array does not match F along dimension $dim. Got $(length(spacing)), must be $(size(F, dim))"
    
    try
        for index ∈ axes(quotient, dim)
            data = selectdim(quotient, dim, index)
            data ./= spacing[index]
        end
        return quotient
    catch e
        @error "Could not compute gradient. $e"
        return nothing
    end
end

const ∇ = gradient