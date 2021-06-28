Data and Snapshot Extraction
============================

In order to get out a snapshot object, simply call the ``snapshot`` function, which is part of the ``dispatch`` module. The arguments are::

        * iout: Integer, snapshot ID
        * run: String, path to snapshots, relative to the data-folder
        * data: String, path to data folder
        * verbose: Int

Example of usage::

        using .dispatch
        s = snapshot(8, data="../data/mhd632")

The resulting object is a dictionairy, as opposed to the class object returned by the Python Implementation. As a result, its
properties (fields) can be accessed using string keys. Examples:: 

        patch1 = s["patches"][1]
        system = s["units"]["system"]
        size = s["cartesian"]["size"]

To get all the available keys in the snapshot, you can do ``println(keys(s))``, or, for prettier output::

        for key in keys(s) println(key) end

Similarly, to iterate through the values, simply call ``values(s)``. To iterate through all key-value pairs, do::

        for (key, val) in s
                do something
        end

A patch is also a dictionairy object, and its properties are accessed in the exact same way. To get 2d-data of a given variable "iv" from a patch, call the ``plane`` function, which takes the arguments::

        * patch: Dict, a patch from a snapshot
        * iv: Int/String, variable name/offset. Default 0
        * x, y, z: Float: coordinate at which to slice the data. Default x=nothing, y=nothing, z=0.5.
        * all: Bool, whether to include guard zones. Default false.

This will produce a 2d-array of type Float32. 

In order to get a buffer of all the data in a snapshot, there are four main functions::

        *``unigrid_plane``
        *``unigrid_volume``
        *``amr_plane``
        *``amr_volume``

If the data you want to extract has been produced using mesh refinent, you have to use the amr-functions, otherwise use unigrid. For the former, you also have to specify an output dimension, i.e. number of cells is each dimension, since the buffer is
is produced by sampling and interpolating. 
