# The `JuliaDispatch.Dispatch` Module

This module contains the `snapshot`-function, which is used for reading a snapshot and parsing the patches.


After having imported the `JuliaDispatch` module into your workspace, simply do

```@repl
using JuliaDispatch.Dispatch
```

The function `snapshot` is then available. See end of page for documentation and examples.

When reading in a snapshot with the `snapshot` function, a dictionary will be returned. All the properties of the snapshot can be accessed using square brackets and a string with the name of the property. For example, if you want to know the time at which the snapshot was taken, do

```julia
println(snap["time"])
```

All the patches and their metadata are stored in the `"patches"` key of the snapshot. This returns a vector with each entry containing the patch metadata in the form of another dictionary. For example

```julia
patch1 = snap["patches"][1]
println(patch1["id"])
```

extracts the first patch and prints its `id`.

The `snapshot` function will by default display a progressbar visualizing the progress in parsing the patches. To disable this, set the keyword argument `progress` to `false`. To suppress all printed information, which might be useful for when reading in a large number of snapshots at once, set `verbose` to `-1`. 

```@docs
JuliaDispatch.Dispatch.snapshot
```


```@docs
JuliaDispatch.Dispatch.cache_snapshots_live
```


```@docs
JuliaDispatch.Dispatch.get_snapshots
```