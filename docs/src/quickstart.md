# Quick-Start Guide


## Installation and Setup

!!! note "Multi-threading"

    `JuliaDispatch` uses multi-threading to speed up the process of parsing patches. It is recommend you start your Julia session with more than one thread. This can be done by typing `julia -threads N` when starting Julia, where `N` is the number of threads. You can also do this using the environment variable by setting `export JULIA_NUM_THREADS=N`. See [docs]https://docs.julialang.org/en/v1/manual/multi-threading/) for more information.

Open your Julia REPL terminal and type ']' to enter the package manager, followed by

```julia
pkg> add "https://github.com/casparwb/JuliaDispatch.jl"
pkg> add JLD
```

`JLD` is a package used by `JuliaDispatch` to cache and read cached namelists. This needs to be in your workspace in order for it work, thus it is important that it is imported before `JuliaDispatch`.

After having installed the packages, exit the package manager with backspace, and type

```julia
using JLD, JuliaDispatch
```

Note that things will be quite slow the first time you add and import a package, in addition to when you call a function. This is due to Julia's precompilation routine, and only happens the first time you call a function.

Importing `JuliaDispatch` does not include any functions into your workspace, but rather allows you to load the individual modules. These modules are

* __Dispatch__: this submodule contains the `snapshot(...)` function which allows reading in a snapshot and saving all the metadata into a dictionary.
* __Graphics__: high-level plotting of slices (planes), volumes, and 1-d quantities.
* __Buffers__: methods for stitching together patch data into 2- or 3-dimensional arrays, in addition to resampling methods for re-sizing the domain.
* __Select__: methods for extracting patches and data at given positions in the computational domain.
* __Utils__:  various QOL tools, such as functions for getting the number of snapshots in a folder, among other things.

After having imported `JuliaDispatch`, you can include the functionality of any of these modules by doing

```julia
using JuliaDispatch.MODULENAME
```

where `MODULENAME` is replaced by any of the abovementioned modules. Note that the module name is case-sensitive, and all modules have a capital first letter. You can also directly use a function within a module without explicitly importing it by doing

```julia
JuliaDispatch.MODULENAME.funcname(args...)
```

where `funcname` is the name of the function exported by the module `MODULENAME`. 

To quickly get the documentation for any function, use the _help_ functionality of the REPL, which you can enter by typing `?`. You can then enter a function for which you want the documentation. For example

```julia
?
snapshot
```

## Functionality

### Snapshots

The `snapshot` function returns a dictionary object. To access any of the snapshot's properties, use square brackets and strings. For example

```julia
time = snapshot["time"]
println(time)
```

gives you the time at which the snapshot was taken. To get all the properties of the snapshot, you can print all the keys or iterate over all key-value pairs

```julia
println(keys(snap))
# or
for (key, value) in snap
    println(key, " ", value)
end
```

### Quantities

To see which quantities are present in the snapshot, you can print the `["idx"]["dict"]` of a snapshot. This returns a dictionary with keys equal to the different quantites, and values as integers corresponding to their offset in the data files. Any quantity that has a value `< 0` is present and can be extracted from the data files.

```@example 1
using JuliaDispatch, JuliaDispatch.Dispatch
data = "../../test/data/orz/data/" # hide
snap = snapshot(0, data=data)
for (k, v) in snap["idx"]["dict"]
    println(k, " ", v)
end
```

### Expression Parsing

The `JuliaDispatch` package has support for parsing expressions. This means that when extracting patch data, instead of requesting a quantity in a specific data slot, a string containing an expression can be sent in instead. The expression parser will then attempt to parse and evaluate the expression and the return the result for the given patch. Any function that accepts the `iv` keyword argument supports this functionality. Variables in your workspace can be interpolated into the expression using the `$` construct. Example

```@example 1
pi = 3.14 #
patch = snap["patches"][1] # get the first patch
velocity = snap["var"]("$pi*sqrt(ux^2 + uy^2 + uz^2)") # equivalent to "3.14*sqrt(ux^2 + uy^2 + uz^2)"
```

### Data Interpolation

One of the main features of the package is the ability to stitch together and interpolate patch data. This can be done with both unigrid and mesh-refined simulations. Additionally, data can be up or downscaled to any given dimensions in order to either increase resolution or save memory. This functionality is exported by the `JuliaDispatch.Buffers` module. See [REF BUFFERS] for details.

## Usage Examples

### Loading a Snapshot

```@example 1
using JuliaDispatch
using JuliaDispatch.Dispatch
data = "../test/data/orz/data" # hide
snap = snapshot(100, data=data) # data is a string variable pointing to the folder containing the snapshots
println("Snapshot time: $(snap["time"])) # get the time the snapshot was taken
println("Number of patches = $(length(snap["patches"]))) # print the number of patches
```

### Extracting Interpolated Data

```@example 1
using JuliaDispatch.Buffers

density_plane = unigrid_plane(snap, iv="d", z=0.1) # get the density in the xy-plane at z=0.1
density_plane_downscaled = amr_plane(snap, iv="d", z=0.1, dims=(200, 300)) # amr_plane can be used to downscale/upscale data

density_volume = unigrid_volum(snap, iv="d") # similar for a 3D volume
density_volume_downscaled = amr_volume(snap, iv="d", dims=100)
```

### Plotting

```@example 1
using JuliaDispatch.Graphics

sliceplot(snap, iv="d", z=1.0) # plot a simple sliceplot of density at z=1.0
sliceplot(snap, iv="sqrt(bx^2+by^2+bz^2)", y=-10) # plot an expression
sliceplot(snap, iv="d", z=0.1, center=(3, 2), width=(6, 6)) # zoom in with (3,2) as the center and a width of (6, 6)
sliceplot(snap, iv="d", z=0.1, resample=true, dims=(400, 600)) # up/down-scale to a size of (400, 600)
sliceplot(snap, iv="ekin", z=0.1, linetype=:surf) # surface plot

plot_values_along(snap, [0.5, 0.5, 0.1], dir=3, iv="d") # plot density along z-axis through point (.5, .5, 0.1)
histogram_along(snap, [0.5, 0.5, 0.1], dir=3, iv="d", norm=:pdf, label="density") # histogram normalized to probability density
anim_pane(snap, ax=3, iv="d", nframes=20, start=0.0, stop=-10.0, savepath="test.gif") # animate a panning in z-direction from -10 to 0 
```