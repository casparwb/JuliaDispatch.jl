
using Unitful

    struct cgs
        name   ::String
        m_earth::Quantity{Float64}
        m_sun  ::Quantity{Float64}
        r_earth::Quantity{Float64}
        grav   ::Quantity{Float64}
        yr     ::Quantity{Float64}
        au     ::Quantity{Float64}
        kms    ::Quantity{Float64}
        mu     ::Quantity{Float64}
        m_u    ::Quantity{Float64}
        k_b    ::Quantity{Float64}
        k_B    ::Quantity{Float64}
        h_p    ::Quantity{Float64}
        h_P    ::Quantity{Float64}
        e      ::Quantity{Float64}
        c      ::Quantity{Float64}
        stefan ::Quantity{Float64}
        m_e    ::Quantity{Float64}
        pc     ::Quantity{Float64}
    end

    struct SI
        name   ::String
        m_earth::Quantity{Float64}
        m_sun  ::Quantity{Float64}
        r_earth::Quantity{Float64}
        grav   ::Quantity{Float64}
        yr     ::Quantity{Float64}
        au     ::Quantity{Float64}
        kms    ::Quantity{Float64}
        mu     ::Quantity{Float64}
        m_u    ::Quantity{Float64}
        m_e    ::Quantity{Float64}
        k_b    ::Quantity{Float64}
        k_B    ::Quantity{Float64}
        h_p    ::Quantity{Float64}
        e      ::Quantity{Float64}
        c      ::Quantity{Float64}
        stefan ::Quantity{Float64}
        pc     ::Quantity{Float64}
    end

    struct scaled
        l::Quantity{Float64}
        d::Quantity{Float64}
        v::Quantity{Float64}
        t::Quantity{Float64}
        m::Quantity{Float64}
        p::Quantity{Float64}
        g::Quantity{Float64}
        u::Quantity{Float64}
        e::Quantity{Float64}
        temp::Quantity{Float64}
    end

    function init_system(units="cgs")
        """ Initialize an object given the system. Returns a struct
        with fields equal to the values """

        if units == "cgs"
            vals = [
            5.972e27u"g",                    # earth mass [g]
            1.989e33u"g",                    # solar mass [g]
            6.371e8u"cm",                    # earth radius [cm]
            6.673e-8u"cm^3 /g / s^2",        # gravitational acceleration [cm³ g⁻¹ s⁻²] 
            3.156e+7u"s",                    # year [s]
            1.498e+13u"cm",                  # AU [cm]
            1e+5u"cm/s",                     # km/s [cm/s]
            2.0u"1.0",                       # mu ???
            1.6726e-24u"g",                  # proton mass [g] 
            1.3807e-16u"erg/K",              # Boltzmann constant [erg K⁻¹]
            1.3807e-16u"g*cm^2/s^2/K",              # Boltzmann constant [erg K⁻¹]
            6.6260e-27u"g*cm^2/s^2/s",              # Planck constant [erg s⁻¹]
            6.6260e-27u"g*cm^2/s^2/s",              # Planck constant [erg s⁻¹]
            4.8032e-10u"sqrt(dyn)cm",              # elementary charge [Statcoulomb]
            2.9979e10u"cm/s",                # speed of light 
            5.6704e-8u"erg / cm^2 /s / K^4", # Stefan-Boltzmann constant [erg cm⁻² s⁻¹ K⁻⁴]
            9.109e-33u"g",                   # electron mass [g]
            3e18u"cm"                        # parsec [cm]
            ]

            return cgs(units, vals...)

        elseif lowercase(units) == "si"
            vals =[
                5.972e24u"kg", # earth mass [kg]
                1.989e30u"kg", # solar mass [kg]
                6.371e6u"m",   # earth radius [m]
                6.673e-11u"m^3 / kg / s^2", # gravitational constant
                3.156e+7u"s", # year [s]
                1.498e+11u"m", # AU [m]
                1e+3u"m/s", # km/s [m/s]
                2.0u"1.0", # mu
                1.6726e-27u"kg", # proton mass
                9.109e-31u"kg", # electron mass
                1.3807e-23u"J/K", # Boltzmann constant [J/K]
                1.3807e-23u"J/K", #  --------=------------
                6.62606e-34u"J / Hz", # Planck constant
                1.6022e-19u"C", # elementary charge [C]
                2.9979e8u"m/s",  # speed of light
                5.6704e-5u"W / m^2 / K^4", # Stefan-Boltzmann constant
                3e16u"m"] # parsec [m]
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
            if system == "cgs"
                l = 1e8u"cm"
                t = 1e2u"s"
                d = 1e-7u"g/cm^3"
                v = (l/t)
            else
                l = 1e8u"m"
                t = 1e2u"s"
                d = 1e-7u"kg/m^3"
                v = (l/t)
            end
        elseif type == "unity"
            l = 1.0u"1.0"
            t = 1e2u"s"
            d = 1.0u"1.0"
            v = l/t
        end
        m = d*l^3
        p = d*v^2
        g = units.grav*d*t^2
        u = units.kms
        e = m*u^2
        temp = units.mu*(units.m_u)/(units.k_b)*v^2

        u_kr = 1.0/(d*l)     # Rosseland opacity
        u_ee = v^2           # specific energy
        u_e = d*u_ee 
        u_te = u_e*v         # box therm. em.
        u_B = v*sqrt(4π*d)
        return scaled(l, d, v, t, m, p, g, u, e, temp)
    end
