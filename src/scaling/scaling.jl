struct cgs
    name   ::String
    m_earth::Float64
    m_sun  ::Float64
    r_earth::Float64
    grav   ::Float64
    yr     ::Float64
    au     ::Float64
    kms    ::Float64
    mu     ::Float64
    m_u    ::Float64
    k_b    ::Float64
    k_B    ::Float64
    h_p    ::Float64
    h_P    ::Float64
    e      ::Float64
    c      ::Float64
    stefan ::Float64
    m_e    ::Float64
    pc     ::Float64
    AU     ::Float64
end

struct SI
    name   ::String
    m_earth::Float64
    m_sun  ::Float64
    r_earth::Float64
    grav   ::Float64
    yr     ::Float64
    au     ::Float64
    kms    ::Float64
    mu     ::Float64
    m_u    ::Float64
    m_e    ::Float64
    k_b    ::Float64
    k_B    ::Float64
    h_p    ::Float64
    e      ::Float64
    c      ::Float64
    stefan ::Float64
    pc     ::Float64
    AU     ::Float64
end

struct scaled
    l::Float64
    d::Float64
    v::Float64
    t::Float64
    m::Float64
    p::Float64
    g::Float64
    u::Float64
    e::Float64
    temp::Float64
end

function init_system(units="cgs")
    """ Initialize an object given the system. Returns a struct
    with fields equal to the values """

    if units == "cgs"
        vals = [
        5.972e27,
        1.989e33,
        6.371e8,
        6.673e-8,
        3.156e+7,
        1.498e+13,
        1e+5,
        2.,
        1.6726e-24,
        1.3807e-16,
        1.3807e-16,
        6.6260e-27,
        6.6260e-27,
        4.8032e-10,
        2.9979e10,
        5.6704e-8,
        9.109e-31,
        3e18,
        1.5e13
        ]

        return cgs(units, vals...)

    elseif lowercase(units) == "si"
        vals =[
            5.972e24,
            1.989e30,
            6.371e6,
            6.673e-11,
            3.156e+7,
            1.498e+11,
            1e+3,
            2.,
            1.6726e-27,
            9.109e-31,
            1.3807e-23,
            1.3807e-23,
            6.62606e-34,
            1.6022e-19,
            2.9979e8,
            5.6704e-5,
            3e16,
            1.5e11]
        return cgs(units, vals...)

    end
end

function scaling(;type="ISM", system="cgs", verbose=0, mu=2)
    """ Returns a struct with scaling constants given the system """

    units = init_system(system)

    verbose > 0 && println("using $(system) units")

    if type == "ISM"
        l = units.pc
        d = 1e-24
        v = 1e5
        t = l/v
    elseif type == "solar"
        l = 1e8
        t = 1e2
        d = 1e-7
        v = l/t
    end

    m = d*l^3
    p = d*v^2
    g = units.grav*d*t^2
    u = units.kms
    e = m*u^2
    temp = mu*(units.m_u)/(units.k_b)*v^2

    return scaled(l, d, v, t, m, p, g, u, e, temp)
end
