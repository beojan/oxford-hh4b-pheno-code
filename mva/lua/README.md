Torch7 MVA implementation
-------------------------

Here we have a very basic implementation of a neural-network MVA in the Torch7
framework. There are two varieties of MVA, detailed in newtorch.lua and
optimtorch.lua. Both of these implement essentially the same MVA but with two
different internal torch procedures.

It should be noted that these scripts have not been extensively benchmarked.

The `generate_test_data.lua` script will generate an example NTuple in the same
format as the normal analysis code. To run one of the MVAs over this test NTuple
just run:

th newtorch.lua test_data.dat
    or
th optimtorch.lua test_data.dat

These will output data files in the same format as the old MVA, which can be
plotted with the plotMVA script found in the ../res folder.


** Better documentation to come **
