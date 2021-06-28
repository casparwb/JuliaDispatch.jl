using ForwardDiff: diff

meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

function curl(axis::Tuple, data::Tuple, dims::Tuple; verbose=0)
    """
    data: tuple of arrays, with field components
    dims: what field components are in data, for example
          data=(Fx, Fy), then dims=(1, 2)
    """

    type = typeof(data[1][1])
    F = Dict{Int, Union{Int, Array{type}}}(1 => 0, 2 => 0, 3 => 0)
    for (f, dim) in zip(data, dims)
        F[dim] = f
    end

    _curl(F[1], F[2], F[3])
end

function _curl(Fx, Fy, Fz)

    dFzdy = 0
    dFydz = 0
    dFxdz = 0
    dFzdx = 0
    dFydx = 0
    dFxdy = 0

    try
        dFzdy = diff(Fz, dims=2)
    catch
        nothing
    end

    try
        dFydz = diff(Fy, dims=3)
    catch
        nothing
    end

    try
        dFxdz = diff(Fx, dims=3)
    catch
        nothing
    end

    try
        dFzdx = diff(Fz, dims=1)
    catch
        nothing
    end

    try
        dFydx = diff(Fy, dims=1)
    catch
        nothing
    end

    try
        dFxdy = diff(Fx, dims=2)
    catch
        nothing
    end

    [dFzdy .- dFydz',
    dFxdz .- dFzdx',
    dFydx .- dFxdy']
end

function divergence(data::Tuple, dims::Tuple; verbose=0)
    """
    data: tuple of arrays, with field components
    dims: what field components are in data, for example
          data=(Fx, Fy), then dims=(1, 2)
    """

    datatype = typeof(data[1][1])
    F = Dict{Int, Union{Int, Array{datatype}}}(1 => 0, 2 => 0, 3 => 0)
    for (f, dim) in zip(data, dims)
        F[dim] = f
    end

    dFxdx = 0
    dFydy = 0
    dFzdz = 0

    try
        dFxdx = diff(F[1], dims=1)
    catch
        nothing
    end

    try
        dFydy = diff(F[2], dims=2)
    catch
        nothing
    end

    try
        dFzdz = diff(F[3], dims=3)
    catch
        nothing
    end

    div = dFxdx .+ dFydy' .+ dFzdz

    return div
end
