""" EOS table reader"""

mutable struct stagger_eos
    md
    iupdate
    nvar
    mbox
    eosxmin
    dbox
    ul
    ut
    ur
    eps
    tff
    grv
    abnd
    tmean
    tamp
    rhm
    xcross
    thmin
    thmax
    dth
    eemin
    eemax
    deetab
    itab
    mtab
    mtable
    tab
    njon
    rhm1
    rhm2
    drhm
end

function initialize(data="table.dat")
    ffreader = pyimport("scipy.io").FortranFile
    ff = ffreader(data, "r")
    (md, iupdate, nvar, mbox, eosxmin),
    (dbox, ul, ut, ur, eps, tff, grv, abnd) = ff.read_record("(5)i4, (8)<f4")[1]

    temp = ff.read_record("(10, $md)f4, (2,$md)i4")[1]
    tmean, tamp, rhm, xcross, thmin, thmax, dth, eemin, eemax, deetab = collect(eachrow(temp[1]))
    itab, mtab = collect(eachrow(temp[2]))

    mtable = ff.read_ints()[1]
    tab = ff.read_record("i4")
    ff.close()

    njon = nvar - 2.0 * mbox
    rhm1 = log(rhm[1])
    rhm2 = log(rhm[end])
    drhm = (rhm2 - rhm1)/(md - 1)

    if njon == 0
        njon = nvar - mbox
    end

    stagger_eos(md, iupdate, nvar, mbox, eosxmin, dbox, ul, ut, ur, eps, tff,
                grv, abnd, tmean, tamp, rhm, xcross, thmin, thmax, dth, eemin,
                eemax, deetab, itab, mtab,mtable, tab, njon, rhm1, rhm2, drhm)

end

function lookup(eos; dim=[20, 20, 20],
                e = nothing, ee = nothing, d = nothing, lnd = nothing,
                pg = nothing, dpg1 = nothing, dpg2 = nothing, tt = nothing,
                ne = nothing, rk = nothing, scr = nothing, ss = nothing)

    mx, my, mz = dim[1], dim[2], dim[3]

    if typeof(lnd) <: AbstractArray
        lnd_loc = lnd
    elseif typeof(d) <: AbstractArray
        lnd_loc = log.(d)
    else
        throw(ArgumentError("no density given to eos!"))
    end

    if typeof(ee) <: AbstractArray
        ee_loc = ee
    elseif e != nothing
        if typeof(lnd) <: AbstractArray
            ee_loc = e ./ d
        else
            ee_loc = e ./ exp.(lnd_loc)
        end
    else
        throw(ArgumentError("no energy given to eos!"))
    end

    for iz ∈ 1:mz
        for iy ∈ 1:my
            nx = zeros(Int, mx)
            px = zeros(mx)

            for ix ∈ 1:mx
                algrk = (lnd_loc[ix, iy, iz] - eos.rhm1)/eos.drhm
                nx[ix] = max(1, min(eos.md-1), Int(round(algrk)))
                px[ix] = algrk - nx[ix]
            end

            ntab = zeros(Int, mx)
            ntab1 = zeros(Int, mx)
            py = zeros(mx)
            ik = zeros(Int, mx)
            ik1 = zeros(Int, mx)

            for ix ∈ 1:mx
                ntab[ix] = eos.mtab[nx[ix]]
                eek = 1.0 + (ee_loc[ix, iy, iz] - eos.eemin[nx[ix]])/eos.deetab[nx[ix]]
                eek = max(eek, 1.0)
                kee = min(ntab[ix]-1, max(1, Int(round(eek))))
                py[ix] = eek - kee
                ik[ix] = eos.itab[nx[ix]] + kee - 1
                ntab1[ix] = eos.mtab[nx[ix]+1]
                kee1 = kee .- Int.(round.(eos.eemin[nx[ix]+1] -
                                          eos.eemin[nx[ix]]) /
                                          eos.deetab[nx[ix]])
                kee1 = min(ntab1[ix]-1, max(1, kee1))
                ik1[ix] = eos.itab[nx[ix]+1] + kee1 - 1
            end

            f00  = zeros(mx)
            f01  = zeros(mx)
            f10  = zeros(mx)
            f11  = zeros(mx)
            fx00 = zeros(mx)
            fx01 = zeros(mx)
            fx10 = zeros(mx)
            fx11 = zeros(mx)
            fy00 = zeros(mx)
            fy01 = zeros(mx)
            fy10 = zeros(mx)
            fy11 = zeros(mx)

            for ix ∈ 1:mx
                qx   = 1.0 - px[ix]
                pxqx = px[ix] * qx
                pxpx = px[ix] * px[ix]
                qxqx = qx * qx

                qy   = 1. - py[ix]
                pyqy = py[ix] * qy
                pypy = py[ix] * py[ix]
                qyqy = qy * qy

                pxqy = px[ix] * qy
                pxpy = px[ix] * py[ix]
                qxqy = qx * qy
                qxpy = qx * py[ix]

                f00[ix] = qxqy * (1. + pxqx - pxpx + pyqy - pypy)
                f01[ix] = qxpy * (1. + pxqx - pxpx + pyqy - qyqy)
                f10[ix] = pxqy * (1. - qxqx + pxqx + pyqy - pypy)
                f11[ix] = pxpy * (1. - qxqx + pxqx + pyqy - qyqy)

                fx00[ix] =    qxqy * pxqx
                fx01[ix] =    qxpy * pxqx
                fx10[ix] =  - pxqy * pxqx
                fx11[ix] =  - pxpy * pxqx

                fy00[ix] =    qxqy * pyqy
                fy01[ix] =  - qxpy * pyqy
                fy10[ix] =    pxqy * pyqy
                fy11[ix] =  - pxpy * pyqy
            end

            if typeof(pg) <: AbstractArray
                println("working on pg")

                for ix ∈ 1:mx
                    pg[ix,iy,iz] = exp(
                    f00[ix] * eos.tab[ik[ix]] +  f01[ix] * eos.tab[ik[ix]  + 1 ]
                 +  f10[ix] * eos.tab[ik1[ix]]+  f11[ix] * eos.tab[ik1[ix] + 1 ]
                 + fx00[ix] * eos.tab[ik[ix]  + 2 * ntab[ix]     ]
                 + fx01[ix] * eos.tab[ik[ix]  + 2 * ntab[ix]  + 1]
                 + fx10[ix] * eos.tab[ik1[ix] + 2 * ntab1[ix]    ]
                 + fx11[ix] * eos.tab[ik1[ix] + 2 * ntab1[ix] + 1]
                 + fy00[ix] * eos.tab[ik[ix]  +     ntab[ix]     ]
                 + fy01[ix] * eos.tab[ik[ix]  +     ntab[ix]  + 1]
                 + fy10[ix] * eos.tab[ik1[ix] +     ntab1[ix]    ]
                 + fy11[ix] * eos.tab[ik1[ix] +     ntab1[ix] + 1]
                                     )
                end # end for
            end # end if

            ik += ntab
            ik1 += ntab1

            if typeof(dpg1) <: AbstractArray
                println("working on dpg1")
                for ix ∈ 1:mx
                    dpg1[ix,iy,iz] = exp(
                    f00[ix] * eos.tab[ik[ix]] +  f01[ix] * eos.tab[ik[ix]  + 1 ]
                 +  f10[ix] * eos.tab[ik1[ix]]+  f11[ix] * eos.tab[ik1[ix] + 1 ]
                 + fx00[ix] * eos.tab[ik[ix]  + 2 * ntab[ix]     ]
                 + fx01[ix] * eos.tab[ik[ix]  + 2 * ntab[ix]  + 1]
                 + fx10[ix] * eos.tab[ik1[ix] + 2 * ntab1[ix]    ]
                 + fx11[ix] * eos.tab[ik1[ix] + 2 * ntab1[ix] + 1]
                 + fy00[ix] * eos.tab[ik[ix]  +     ntab[ix]     ]
                 + fy01[ix] * eos.tab[ik[ix]  +     ntab[ix]  + 1]
                 + fy10[ix] * eos.tab[ik1[ix] +     ntab1[ix]    ]
                 + fy11[ix] * eos.tab[ik1[ix] +     ntab1[ix] + 1]
                                       )
                end # end for
            end # end if

            ik += ntab
            ik1 += ntab1

            if typeof(rk) <: AbstractArray
                println("working on rk")
                for ix ∈ 1:mx
                    rk[ix,iy,iz] = exp(
                    f00[ix] * eos.tab[ik[ix]] +  f01[ix] * eos.tab[ik[ix]  + 1 ]
                 +  f10[ix] * eos.tab[ik1[ix]]+  f11[ix] * eos.tab[ik1[ix] + 1 ]
                 + fx00[ix] * eos.tab[ik[ix]  + 2 * ntab[ix]     ]
                 + fx01[ix] * eos.tab[ik[ix]  + 2 * ntab[ix]  + 1]
                 + fx10[ix] * eos.tab[ik1[ix] + 2 * ntab1[ix]    ]
                 + fx11[ix] * eos.tab[ik1[ix] + 2 * ntab1[ix] + 1]
                 + fy00[ix] * eos.tab[ik[ix]  +     ntab[ix]     ]
                 + fy01[ix] * eos.tab[ik[ix]  +     ntab[ix]  + 1]
                 + fy10[ix] * eos.tab[ik1[ix] +     ntab1[ix]    ]
                 + fy11[ix] * eos.tab[ik1[ix] +     ntab1[ix] + 1]
                                     )
                end # end for

                if size(rk)[4] > 1
                    for ij ∈ 2:size(rk)[4]
                        rk[:, iy, iz, ij] = rk[:, iy, iz, ij-1]*10
                    end # end for
                end
            end # end if

            ik += 3*ntab
            ik1 += 3*ntab1

            if typeof(tt) <: AbstractArray
                println("working on tt")
                for ix ∈ 1:mx
                    tt[ix,iy,iz] = exp(
                    f00[ix] * eos.tab[ik[ix]] +  f01[ix] * eos.tab[ik[ix]  + 1 ]
                 +  f10[ix] * eos.tab[ik1[ix]]+  f11[ix] * eos.tab[ik1[ix] + 1 ]
                 + fx00[ix] * eos.tab[ik[ix]  + 2 * ntab[ix]     ]
                 + fx01[ix] * eos.tab[ik[ix]  + 2 * ntab[ix]  + 1]
                 + fx10[ix] * eos.tab[ik1[ix] + 2 * ntab1[ix]    ]
                 + fx11[ix] * eos.tab[ik1[ix] + 2 * ntab1[ix] + 1]
                 + fy00[ix] * eos.tab[ik[ix]  +     ntab[ix]     ]
                 + fy01[ix] * eos.tab[ik[ix]  +     ntab[ix]  + 1]
                 + fy10[ix] * eos.tab[ik1[ix] +     ntab1[ix]    ]
                 + fy11[ix] * eos.tab[ik1[ix] +     ntab1[ix] + 1]
                                      )
                end # end for
            end # end if

            ik += 3*ntab
            ik1 += 3*ntab1

            if typeof(ne) <: AbstractArray
                println("working on ne")
                for ix ∈ 1:mx
                    ne[ix,iy,iz] = exp(
                     f00[ix] * eos.tab[ik[ix]] +  f01[ix] * eos.tab[ik[ix]  + 1 ]
                  +  f10[ix] * eos.tab[ik1[ix]]+  f11[ix] * eos.tab[ik1[ix] + 1 ]
                  + fx00[ix] * eos.tab[ik[ix]  + 2 * ntab[ix]     ]
                  + fx01[ix] * eos.tab[ik[ix]  + 2 * ntab[ix]  + 1]
                  + fx10[ix] * eos.tab[ik1[ix] + 2 * ntab1[ix]    ]
                  + fx11[ix] * eos.tab[ik1[ix] + 2 * ntab1[ix] + 1]
                  + fy00[ix] * eos.tab[ik[ix]  +     ntab[ix]     ]
                  + fy01[ix] * eos.tab[ik[ix]  +     ntab[ix]  + 1]
                  + fy10[ix] * eos.tab[ik1[ix] +     ntab1[ix]    ]
                  + fy11[ix] * eos.tab[ik1[ix] +     ntab1[ix] + 1]
                                      )
                end # end for
            end # end if

            if typeof(scr) <: AbstractArray
                println("working on scr")
                for ix ∈ 1:mx
                    scr[ix,iy,iz] = exp(
                    f00[ix] * eos.tab[ik[ix]] +  f01[ix] * eos.tab[ik[ix]  + 1 ]
                 +  f10[ix] * eos.tab[ik1[ix]]+  f11[ix] * eos.tab[ik1[ix] + 1 ]
                 + fx00[ix] * eos.tab[ik[ix]  + 2 * ntab[ix]     ]
                 + fx01[ix] * eos.tab[ik[ix]  + 2 * ntab[ix]  + 1]
                 + fx10[ix] * eos.tab[ik1[ix] + 2 * ntab1[ix]    ]
                 + fx11[ix] * eos.tab[ik1[ix] + 2 * ntab1[ix] + 1]
                 + fy00[ix] * eos.tab[ik[ix]  +     ntab[ix]     ]
                 + fy01[ix] * eos.tab[ik[ix]  +     ntab[ix]  + 1]
                 + fy10[ix] * eos.tab[ik1[ix] +     ntab1[ix]    ]
                 + fy11[ix] * eos.tab[ik1[ix] +     ntab1[ix] + 1]
                                      )
                end # end for
            end # end if
        end # end iy
    end # end iz
end # end lookup function
