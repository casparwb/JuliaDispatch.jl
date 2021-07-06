<!-- # Overview

_JuliaDispatch.jl_ is a package developed for analysis and visualization of simulations from the Dispatch framework.

## Installation

To install the package, enter the Julia package manager (Pkg) by typing `]` in the Julia REPL, followed by

`pkg> add "https://github.com/casparwb/JuliaDispatch.jl"`

You can check that your environment contains the package by typing (still inside Pkg)

`pkg> status`

which should list `JuliaDispatch` as one of the available packages. To utilize the package, import it into your workspace by typing

`julia> using JuliaDispatch`

This might take a while as Julia has to precompile the package. 

## Getting Started

The _JuliaDispatch_ package contains multiple submodules each of which contains functions and methods for performing various tasks, such as reading in a snapshot, visualization, and data buffering. The full list of modules is

* _Dispatch_: this submodule contains the `snapshot(...)` function which allows reading in a snapshot and saving all the metadata into a dictionary.
* _Graphics_: high-level plotting of slices (planes), volumes, and 1-d quantities.
* _Buffers_: methods for stitching together patch data into 2- or 3-dimensional arrays, in addition to resampling methods for re-sizing the domain.
* _Select_: methods for extracting patches and data at given positions in the computational domain.
* DispatchUtils:  -->