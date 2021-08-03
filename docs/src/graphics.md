# The `JuliaDispatch.Graphics` Module

This module contains various functions for visualizations, including 1D, 2D and 3D plots and simple animations.

To begin, first import the module into your workspace 

```@repl
using JuliaDispatch.Graphics
```

Assuming you already have a snapshot loaded, you can then call the various plotting functions to visualize quantities. All the plotting functions accepts keyword arguments which are supported by the `Plots.jl` GR backend. See [documentation](http://docs.juliaplots.org/latest/generated/gr/) for all possible keyword arguments. This means that you can set title, axis labels, fontsizes, linestyle (lineplot, scatter, heatmap, surface, etc) manually. 

## `sliceplot`

The `sliceplot` function is the main function for plotting a quantity in a given plane. The function is very flexible and includes, among other things, the possibility to zoom in on certain regions by using the `center` and/or `width` keyword arguments. If the given combination of center and width extends beyond a boundary, the function will look for periodicity in the snapshot and stitch together the data beyond the boundary with periodic axes.

By default, `sliceplot` assumes the snapshot data is not mesh-refined. If it is mesh-refined, make sure the `unigrid` keyword argument is set to `false`.

It it possible to rescale the data by setting the `dims` keyword argument to either an integer (which will scale the domain to `NxN`), or a length-2 collection, which scales the domain to `NxM`. 

```@docs
JuliaDispatch.Graphics.sliceplot
```

## `plot_values_along`

To plot values of a quantity along a given axis through a given point, use the function `plot_values_along`. This function takes in a snapshot, a point (in the form of a length-3 array) as arguments, and the quantity `iv` and the axis direction `dir` as keyword arguments, in addition to any keyword argument accepted by the `Plots.plot()` function.

```@docs
JuliaDispatch.Graphics.plot_values_along
```

## `histogram_along`

The function `histogram_along` computes and plots a histogram of a quantity `iv` along a direction `dir` through a point `pt`. By default the histogram shows the absolute frequency (counts), but can be changed to a probability density or probability by setting the keyword argument `norm` to, respectively, `:pdf` or `:probability`. Additionally you define the bins with the `bins` keyword argument. 

```@docs
JuliaDispatch.Graphics.histogram_along
```

## `anim_plane`

To produce a time-evolution animation of given quantity in a 2D plane, use the `anim_plane` function. This function accepts any keyword arguments accepted by `sliceplot`, in addition to the `tspan` keyword argument which can be a tuple of numbers denoting at which time (which snapshot) to start and end the simulation. Additionally, the `step` keyword arguments defines whether to use all snapshots or to skip with a given `step` value.

```@docs
JuliaDispatch.Graphics.anim_plane
```

## `anim_pane`

If you want animate a pane through the computational domain in a certain direction, call the `anim_pane` function. This function accepts any keyword argument supported by `sliceplot` in addition to the `nframes` keyword argument which denotes how many frames will be written to the output video file. By default a pan will be made from the minimum value of the given axis to the maximum value, but this can be overwritten by setting the keyword arguments `start` and/or `stop` to desired values. Additionally you can set `reverse` to `true` if you would like to pan from the maximum value to the minimum. 

```@docs
JuliaDispatch.Graphics.anim_pane
```

## `plot_time_evolution`

```@docs
JuliaDispatch.Graphics.plot_time_evolution
```