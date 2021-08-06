# using ForwardDiff: diff

"""
Compute the approximate curl given the components of a vector field.

Input
- Fs: tuple or array of arrays; components of vector field
- spacing: number or array; spacing between values in the components, for
         computing the gradient. See gradient() function for more details.


Output
- 2- or 3-dimensional array of arrays with curl of field in given
  dimensions.
"""
function curl(Fs::Union{Tuple, AbstractArray},
              spacing::Union{Number, AbstractArray}=1)

    if length(Fs) == 3

        F1, F2, F3 = Fs
        dF3d2 = gradient(F3, dims=2, spacing=spacing)
        dF2d3 = gradient(F2, dims=3, spacing=spacing)

        dF1d3 = gradient(F1, dims=3, spacing=spacing)
        dF3d1 = gradient(F3, dims=1, spacing=spacing)

        dF2d1 = gradient(F2, dims=1, spacing=spacing)
        dF1d2 = gradient(F1, dims=2, spacing=spacing)

        return [dF1d2 .- dF2d3,
                dF1d3 .- dF3d1,
                dF2d1 .- dF1d2]
    else
        F1, F2 = Fs

        dF2d1 = gradient(F2, dims=1, spacing=spacing)
        dF1d2 = gradient(F1, dims=2, spacing=spacing)

        return [dF1d2,
                dF2d1 .- dF1d2]
    end

end

"""

Function for computing the divergence of a vector field given components.

Input

- Fs: tuple or array or arrays; components of field for which to find the
    divergence
- spacing: number or array; spacing between values in the components, for
         computing the gradient. See gradient() function for more details.
         Default 1.

Output

- array with divergence of given components
"""
function divergence(Fs::Union{Tuple, AbstractArray},
                    spacing::Union{Number, AbstractArray}=1)

    dims = collect(1:ndims(Fs[1]))
    grads = [gradient(F, dims=dim, spacing=spacing) for (F, dim) in zip(Fs, dims)]

    result = zeros(eltype(grads[1]), size(grads[1])...)
    for i = 1:length(grads)
        result = result .+ grads[i]
    end
    return result
end

"""
Gradient of an N-dimensional array along given dimensions

Input
- F: N-dimensional array of any numeric type; values of which to differentiate
- dims: nothing, Integer, Array or Tuple; which dimension(s) to differentiate
        F along
- spacing: Number or Array; spacing between the values in F. If Number,
           uniform spacing is used on all values in all dimensions. If
           1D-array, must be length N, uniform distance for each dimension.
           If ND-array, spacing for each value along each dimension.

Output
- array or array of arrays with gradients along each given dim in dims
  (1d if dims is Number)
"""
function gradient(F;
                  dims::Union{Nothing, Int, Array, Tuple}=nothing,
                  spacing::Union{Number, AbstractArray}=1)

    if any(dims .> ndims(F))
        println("F does not contain dimension")
        return 0
    end

    grad = _gradient(F, dims=dims, spacing=spacing)

    if length(grad) == 1
        return grad[1]
    else
        return grad
    end
end


function _gradient(F;
                   dims::Union{Nothing, Int, Array, Tuple}=nothing,
                   spacing::Union{Number, Array, Tuple}=1)

    N = nothing
    if isnothing(dims)
        N = 1:ndims(F)
    else
        N = dims
    end

    isa(N, Int) ? N = [N] : nothing
    Ndims = length(N)
    Fsize = size(F)

    dF = [similar(F) for i = 1:Ndims]
    # dF = [diff(F, dims=dim) for dim in N]
    for i = 1:Ndims
        dF_ = diff(F, dims=N[i])
        sze = [size(dF_)...]
        setindex!(sze, 1, N[i])
        dF[i] = cat(dF_, zeros(sze...), dims=N[i])
    end

    if typeof(spacing) <: Number
        return dF/spacing
    elseif typeof(spacing) <: AbstractArray
        if ndims(spacing) > 1
            @assert ndims(spacing) == Ndims
                    "Number of spacing arrays must match input array dimensions"

        return [dFf ./ h for (dFf, h) in zip(dF, spacing)]
        else
            @assert length(spacing) == Ndims
                    "Number of spacing scalars must match input array dimensions"
            return [dFf / h for (dFf, h) in zip(dF, spacing)]
        end
    end
end
