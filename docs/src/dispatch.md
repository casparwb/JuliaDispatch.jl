# The `JuliaDispatch.Dispatch` Module

This module contains the `snapshot`-function, which is used for reading a snapshot and parsing the patches.

## Usage

After having imported the `JuliaDispatch` module into your workspace, simply do

`julia> using JuliaDispatch.Dispatch`

The function `snapshot` is then available. See end of page for documentation and examples.

When reading in a snapshot with the `snapshot` function, a dictionary will be returned. All the properties of the snapshot can be accessed using square brackets and a string with the name of the property. For example, if you want to know the time at which the snapshot was taken, do

`julia> snap["time"]`



## Function Documentation

```@docs
JuliaDispatch.Dispatch.snapshot
```
