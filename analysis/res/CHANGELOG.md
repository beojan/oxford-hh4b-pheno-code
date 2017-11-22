CHANGELOG
==============

Key Changes
-----------

### Merging

- `yodamerge` rewritten in C++ to avoid needing to load Python and `yoda` module
- New `yoda-merge` also writes out a flat file, removing need to run `yoda2flat`
- Merging scripts ported to ZSH. Compiling with `zcompile` speeds these up
- In total, merging is sped up by a factor of 4 to 6

Unfortunately, `yoda-merge.cpp` doesn't currently build on PPLXINT (update: loading LCG 91 may help).
Either use the binary in the git repo, or build it with `gcc -o yoda-merge yoda-merge.cpp -O2 -lboost_program_options`
on an x86_64 (i.e. any recent) Linux system with a recent GCC and Boost. I've only tested with
GCC 6.3 and Boost 1.63, but it *should* work with slightly older versions (likely GCC 5 and above,
though adding `-std=gnu++14` may be advisable on versions older than 6.1).

### Plotting

The main factor causing plotting to be slow was the time taken to load Python, and the Numpy and
Matplotlib modules. This can be mitigated by doing all plotting in one Python script. The `plotting`
module helps with this.

- New `plotting` Python modules has been added providing two functions:
  - `plot1D` for plotting (multiple) 1D histograms, normalized by area.
  - `plot2D` for plotting a single 2D histogram

- New `plotting_matchpaper` Python modules has been added which can be used to produce plots matching the previous paper:
  - It is the same as `plotting` except bin content is not multiplied by bin width.

- New `makeplots_matchpaper` Python modules has been added as an example of using the functions in plotting_matchpaper.


### Environment Setup

The code requires some environment setup. Much of the Python related setup has
been moved into a virtualenv using Python 2.7.13 built with GCC 4.9.3, and with
Numpy, Matplotlib, and Seaborn (needed for the plotting module) installed.

All environment setup can be done on the PPLXINT machines by sourcing the setupEnv.sh script in the package.
