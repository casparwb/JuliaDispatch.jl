"""
   Support functions for reading or calculating DISPATCH2 grid and geometry information.

"""

function GeometricFactors(p)
    ### Constructor ###
    GF = Dict(
            "h2c" => nothing,
            "h2f" => nothing,
            "h31c" => nothing,
            "h31f" => nothing,
            "h32c" => nothing,
            "h32f" => nothing,
            "dx1c" => nothing,
            "dx1f" => nothing,
            "dx2c" => nothing,
            "dx2f" => nothing,
            "dx3c" => nothing,
            "dx3f" => nothing,
            "dvol1c" => nothing,
            "dvol1f" => nothing,
            "dvol2c" => nothing,
            "dvol2f" => nothing,
            "dvol3c" => nothing,
            "dvol3f" => nothing,
            "dar1c" => nothing,
            "dar1f" => nothing,
            "dar2c" => nothing,
            "dar2f" => nothing,
            "dar31c" => nothing,
            "dar31f" => nothing,
            "dar32c" => nothing,
            "dar32f" => nothing
            )

    init_grid(GF, p)

end

function init_grid(GF, p)
    """Initialise geometric factors based on coord. type."""

    if p["mesh_type"] == "Cartesian"
        init_Cartesian(GF, p)
    elseif p["mesh_type"] == "cylindrical"
        init_cylindrical(GF, p)
    elseif p["mesh_type"] == "spherical"
        init_spherical(GF, p)
    end

end

function init_Cartesian(GF, p)
    """Initialise geometric factors for a Cartesian coord. system."""

    n1, n2, n3 = p["ncell"]

    # 1-direction
    GF["h2c"] = ones(n1)
    GF["h2f"] = ones(n1)
    GF["h31c"] = copy(GF["h2c"]) #.view()
    GF["h31f"] = copy(GF["h2f"])#.view()

    # 2-direction
    GF["h32c"] = ones(n2)
    GF["h32f"] = copy(GF["h32c"])#.view()

    # linear size elements
    GF["dx1c"] = ones(n1) * p["ds"][1]
    GF["dx1f"] = ones(n1) * p["ds"][1]
    GF["dx2c"] = ones(n2) * p["ds"][2]
    GF["dx2f"] = ones(n2) * p["ds"][2]
    GF["dx3c"] = ones(n3) * p["ds"][3]
    GF["dx3f"] = ones(n3) * p["ds"][3]

    # volume elements
    GF["dvol1c"] = ones(n1) * p["ds"][1]
    GF["dvol1f"] = ones(n1) * p["ds"][1]
    GF["dvol2c"] = ones(n2) * p["ds"][2]
    GF["dvol2f"] = ones(n2) * p["ds"][2]
    GF["dvol3c"] = ones(n3) * p["ds"][3]
    GF["dvol3f"] = ones(n3) * p["ds"][3]

    # area elements
    GF["dar1c"] = GF["h31c"] * GF["h2c"]
    GF["dar1f"] = GF["h31f"] * GF["h2f"]
    GF["dar2c"] = GF["h31f"] * p["ds"][1] / GF["dvol1c"]
    GF["dar2f"] = GF["h31c"] * p["ds"][1] / GF["dvol1f"]
    GF["dar31c"] = GF["h2f"] * p["ds"][1] / GF["dvol1c"]
    GF["dar31f"] = GF["h2c"] * p["ds"][1] / GF["dvol1f"]
    GF["dar32c"] = p["ds"][2] / GF["dvol2c"]
    GF["dar32f"] = p["ds"][2] / GF["dvol2f"]


end

function init_cylindrical(GF, p)
    """ Initialise geometric factors for a cylindrical coordinate system """

    n1, n2, n2 = p["cell"]

    # 1-direction
    GF["h2c"] = ones(n1)
    GF["h2f"] = ones(n1)
    GF["h31c"] = view(GF["h2c"])
    GF["h31f"] = view(GF["h2f"])

    # 2-direction
    pos_c = copy(p["y"])
    pos_f = copy(p["ys"])
    GF["h32c"] = abs.(pos_c)
    GF["h32f"] = abs.(pos_f)

    # linear size elements
    GF["dx1c"] = ones(n1) .* p.ds[1]
    GF["dx1f"] = ones(n1) .* p.ds[1]
    GF["dx2c"] = ones(n2) .* p.ds[2]
    GF["dx2f"] = ones(n2) .* p.ds[2]
    GF["dx3c"] = ones(n3) .* p.ds[3]
    GF["dx3f"] = ones(n3) .* p.ds[3]

    # volume elements
    GF["dvol1c"] = ones(n1) .* p["ds"][1]
    GF["dvol1f"] = ones(n1) .* p["ds"][1]
    GF["dvol2c"] = similar(pos_c)
    GF["dvol2f"] = similar(pos_f)

    for j = 1:length(pos_c)
        GF["dvol2c"][j] = 0.5 * abs(GF["h32f"][j+1] * pos_f[j+1]
                                  - GF["h32f"][j  ] * pos_f[j  ] )
        if j == 1
            GF["dvol2f"][j] = 0.5 * abs(GF["h32c"][j  ] * pos_c[j]
                                      - GF["h32c"][end] * pos_c[end] )
        else
            GF["dvol2f"][j] = 0.5 * abs(GF["h32c"][j  ] * pos_c[j]
                                      - GF["h32c"][j-1] * pos_c[j-1] )
        end
    end

    GF["dvol3c"] = ones(n3) .* p["ds"][2]
    GF["dvol3f"] = ones(n3) .* p["ds"][2]

    # area elements
    GF["dar1c"]  = GF["h31c"] .* GF["h2c"]
    GF["dar1f"]  = GF["h31f"] .* GF["h2f"]
    GF["dar2c"]  = GF["h31f"] .* p.ds[1] ./ GF["dvol1c"]
    GF["dar2f"]  = GF["h31c"] .* p.ds[1] ./ GF["dvol1f"]
    GF["dar31c"] = GF["h2f"]  .* p.ds[1] ./ GF["dvol1c"]
    GF["dar31f"] = GF["h2c"]  .* p.ds[1] ./ GF["dvol1f"]
    GF["dar32c"] = p.ds[2] ./ GF["dvol2c"]
    GF["dar32f"] = p.ds[2] ./ GF["dvol2f"]

end

function init_spherical(GF, p)
    """ Initialise geometric factors for a spherical coord. system """

    n1, n2, n3 = p["ncell"]

    # 1-direction
    rpos_c = copy(p["x"])
    rpos_f = copy(p["xs"])
    GF["h2c"] = abs.(rpos_c)
    GF["h2f"] = abs.(rpos_f)
    GF["h31c"] = view(GF["h2c"])
    GF["h31f"] = view(GF["h2f"])

    # 2-direction
    tpos_c = copy(p["y"])
    tpos_f = copy(p["ys"])
    GF["h32c"] = abs.(sin.(tpos_c))
    GF["h32f"] = abs.(sin.(tpos_f))

    # linear size elements
    GF["dx1c"] = ones(n1) .* p["ds"][1]
    GF["dx1f"] = ones(n1) .* p["ds"][1]
    GF["dx2c"] = ones(n2) .* p["ds"][2]
    GF["dx2f"] = ones(n2) .* p["ds"][2]
    GF["dx3c"] = ones(n3) .* p["ds"][3]
    GF["dx3f"] = ones(n3) .* p["ds"][3]

    # volume elements
    GF["dvol1c"] = similar(rpos_c)
    GF["dvol1f"] = similar(rpos_f)
    for i = 1:length(rpos_c)

        GF["dvol1c"] = (GF["h2f"][i+1] * GF["h31f"][i+1] * rpos_f[i+1] / 3.0
                      - GF["h2f"][i  ] * GF["h31f"][i  ] * rpos_f[i  ] / 3.0)
        if i == 1
            GF["dvol1f"] = (GF["h2c"][i  ] * GF["h31c"][i  ] * rpos_c[i  ] / 3.0
                         -  GF["h2c"][end] * GF["h31c"][end] * rpos_c[end] / 3.0)
        else
            GF["dvol1f"] = (GF["h2c"][i  ] * GF["h31c"][i  ] * rpos_c[i  ] / 3.0
                         -  GF["h2c"][i-1] * GF["h31c"][i-1] * rpos_c[i-1] / 3.0)
        end
   end

    GF["dvol2c"] = similar(tpos_c)
    GF["dvol2f"] = similar(tpos_f)

    for j = 1:length(tpos_c)
        GF["dvol2c"] = -cos(tpos_c[j+1]) + cos(tpos_c[j  ])
        if i == 1
            GF["dvol2f"] = -cos(tpos_f[j  ]) + cos(tpos_f[end])
        else
            GF["dvol2f"] = -cos(tpos_f[j  ]) + cos(tpos_f[j-1])
        end
    end

    GF["dvol3c"] = ones(n3) .* p.ds[3]
    GF["dvol3f"] = ones(n3) .* p.ds[3]

    # area elements
    GF["dar1c"] = GF["h31c"] .* GF["h2c"]
    GF["dar1f"] = GF["h31f"] .* GF["h2f"]
    GF["dar2c"] = GF["h31f"] .* p.ds[1] ./ GF["dvol1c"]
    GF["dar2f"] = GF["h31c"] .* p.ds[1] ./ GF["dvol1f"]
    GF["dar31c"] = GF["h2f"] .* p.ds[1] ./ GF["dvol1c"]
    GF["dar31f"] = GF["h2c"] .* p.ds[1] ./ GF["dvol1f"]
    GF["dar32c"] = p.ds[2] ./ GF["dvol2c"]
    GF["dar32f"] = p.ds[2] ./ GF["dvol2f"]

end # end init_spherical()
