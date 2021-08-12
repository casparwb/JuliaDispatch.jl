"""
Example usage of JuliaDispatch.jl using data from the Kelvin-Helmholtz Instability
experiment
"""

cd(@__DIR__) # move to current folder

### Setting up the environment ###
begin 
    using Pkg           # import the package manager
    Pkg.activate(".")   # use this folder as environment
    Pkg.add(url="https://github.com/casparwb/JuliaDispatch.jl")
    Pkg.add("Plots")
    Pkg.update()

    using JuliaDispatch, Plots
    using JuliaDispatch.Dispatch, JuliaDispatch.Utils, JuliaDispatch.Graphics, JuliaDispatch.Analysis, JuliaDispatch.Buffers
end

### Get the snapshot to work with ###
begin
    data = "data/"      # directory containing data
    rundir = "perturb"     # name of simulation run

    snapIDs = get_snapshot_ids(data=data, run=rundir)   # snapshot IDs
    
    snap = snapshot(snapIDs[end], data=data, run=rundir); # get the first snapshot
end

### Plotting ###
begin
    p1 = sliceplot(snap, iv="d", z=0, title="""density @ t=$(snap["time"]) @ z=0""", clims=(0.9, 2.1))  # produce a sliceplot of the density at z = 0
    p2 = sliceplot(snap, iv="d", z=0, span=((0, 0.2), (-0.2, 0.0)), clims=(0.9, 2.1)) # plot the density at x ∈ [0.0, 0.2], y ∈ [-0.2, 0.0], z = 0.0
    p3 = plot_values_along(snap, [0.0, 0.1, 0.0], dir=1, iv="d")    # plot the density along the x-axis through (0.0, 0.1, 0.0)

    plot(p1, p2, p3, layout=@layout([a b; c]), size=(800, 800))      # combine all three plots into one subplot
    savefig("density.svg") # save to svg file in current directory
end


### Animating ###
begin
    anim_plane(data=data, run=run, iv="d", z=0, clims=(0.0, 0.5), savepath="density_z=0.gif") # animate time evolution of density at z=0
    anim_plane(data=data, run=run, iv="d", z=0, clims=(0.0, 0.5), tspan=(5, 7), savepath="density_z=0.gif") # animate only from t=5 to t=7

    anim_pane(snap, ax=3, iv="d", start=-20, stop=-4, savepath="density_along_z.gif", clims=(300, 30_000)) # animate density from z=-20 to z=-4
end


### Analysis ###
begin
    average(snap, iv="d") # average of density in whole computational domain
    time_evolution_average(iv="d", data=data, run=run) # time evolution of density average
    plaverage(snap, iv="d", dir=1, nslices=100) # average of density in 100 yz-planes along x-axis 

    dens = get_plane(snap, iv="d", z=0); # 2D array containing all density data

    Δ = snap["ds"][1] # spacing in x- and y-axes

    d_dens_dx, d_dens_dy = ∇(dens, spacing=Δ[1]) # gradient of density 

    dens_gradient == gradient(dens, spacing=Δ[1]) # ∇ is an alias for gradient

    heatmap(d_dens_dx + d_dens_dy)
end