# The `JuliaDispatch.Buffers` Module

The purpose of this module is to _stitch_ together all the patch data of a given quantity in to a single array. The module contains functionality for doing this for both a unigrid and a mesh-refined simulation, and for both a 3D and 2D domain. Additionally, it has the ability to down- or upscale the stitched data using gridded bi/tri-linear interpolation, which allows for faster and more memory efficient data extraction.

As with the other modules, after having imported ´JuliaDispatch` into your workspace, do

```@repl
using JuliaDispatch.Buffers
```

## `unigrid_plane`

The `unigrid_plane` function takes in a snapshot and a single physical coordinate at some `x/y/z`-value and attempts to combine the data using the patches that are found within the given plane. This function returns the entire data set, meaning it does not perform any interpolation when stitching together the data apart from interpolating patch data between planes to get the values at the given coordinate.

```@docs
JuliaDispatch.Buffers.unigrid_plane
```

## `unigrid_volume`

Similar to `unigrid_plane` except that it returns a 3D array with the entire simulation data set for the given quantity.

```@docs
JuliaDispatch.Buffers.unigrid_volume
```

## `amr_plane`

If the simulation uses mesh-refinement __or__ you want to up/downscale the data from a unigrid experiment, use the `amr_plane` function. In addition to the standard parameters of `unigrid_plane` , this function also accepts the dimensions of the returned array containing the interpolated data. The dimensions are given with the keyword argument `dims` and can be an `Integer` denoting a square array, or a length-2 array or tuple denoting the dimensions in axis 1 and axis 2. 

```@docs 
JuliaDispatch.Buffers.amr_plane
```

## `amr_volume`

Similar to `amr_volume` but for 3D data.

```@docs
JuliaDispatch.Buffers.amr_volume
```


