a:

Getting Started
===================
To use the Julia DISPATCH analysis tools, you first need to include all the modules. This can be done by including the ``modules.jl``-file, which is located in ``utilities/julia``. In the REPL terminal, run::

    include(path_to_file)    

For example, if you are located in a folder ``experiments/experimentname/julia``, you would do::

    include("../../../utilities/julia/modules.jl")

After doing this, you can import the individual modules, of which there are 4 available; **dispatch**, **graphics**, **select**, and **analysis**. To import any of them, type::
    
    using .module_name

where module_name is any of the abovementioned. Remember to include the "." in front of the name. You can also do ``import`` instead of ``using``, wherein the difference is that when doing import you need to explicitly prefix the functions with the module the are imported from, whereas when doing ``using`` they are directly added to the namespace. For example::

    import .dispatch
    snap = dispatch.snapshot(8)

or::

    using .dispatch
    snap = snapshot(8)

If you want to use prefixes, but prefer something shorter, you can set your own by doing::

    import .dispatch
    const prefix = dispatch

where prefix is a constant with whatever name you want. 

Exported Functions
==================

Not all functions in the modules are exported, meaning that even if you import them with ``using``, you still need to prefix them. Following is a list of all functions in each module, and which of these are exported.

* dispatch

    * Exported:

        * snapshot


Extracting snapshots and data
=============================


