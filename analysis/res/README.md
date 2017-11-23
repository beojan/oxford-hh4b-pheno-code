README
==============

### Merging

Merging of the produced histogram.yoda files and the ntuple.dat files into histogram.dat files and merged ntuple.dat files is performed on the batch system using the merge.job script. The default is to merge files in directory `res/baseline_noPU_atlas_qcd`, however, the script can be used to merge files in other output directories by changing all instances of `baseline_noPU_atlas_qcd` in the script to the desired directory name (e.g. `baseline_noPU_oxford`). Note that the shebang line in the script is `#!/usr/bin/env zsh` as the `merge.job` script calls the `mergeSubSamples.zsh` script. The `mergeSubSamples.zsh` script utilises yoda-merge (this is pre-compiled, however, if you need to re-compile it please see the details in the CHANGELOG). To run the script on the batch use the following command:

```
qsub merge.job
```

Following this, the mergeBackground.zsh script is run from the output directory (e.g. `baseline_noPU_atlas_qcd`) to combine all histograms and ntuples in the SHERPA_QCD2b2j, SHERPA_QCD4b and SHERPA_QCD4j directories into combined background histograms and ntuples in a new `background` directory (e.g. `baseline_noPU_atlas_qcd/background`). In order to run this script do the following:

```
cd baseline_noPU_atlas_qcd 
```
or directory name of your choosing. 

```
../mergeBackground.zsh
```

If running the `oxford_atlas_qcd` algorithm (e.g. for background estimation studies) an additional merging step is required in order to merge the produced `fullNTuple` dat files, executed as follows:

```
cd baseline_noPU_atlas_qcd 
```

```
../mergeFullNTuple.zsh
```

Note that to speed up merging you can compile the appropriate merging script (`mergeSubSamples.zsh`, `mergeBackground.zsh`  or `mergeFullNTuple.zsh`) by opening a Z shell via `zsh`, compiling with `zcompile mergeSubSamples.zsh` for example, then executing the script as described above.

### Plotting

There are two plotting modules `plotting` and `plotting_matchpaper` which can be imported via `from plotting import *` or `from plotting_matchpaper import *`. Each contains `plot1D` and `plot2D`
functions which can be used. `plotting_matchpaper` is the same as `plotting` except that bin content is not multiplied by bin width, producing plots which match the previous paper. An example of using the functions in `plotting_matchpaper` is given in `makeplots_matchpaper.py`. Note that assuming the environment has been setup as described previously (i.e. you are using the virtualenv), the shebang line for the `plotting_matchpaper` script should be `#!/usr/bin/env python2`. Do *not* use `#!/usr/bin/python` as this will use the system Python interpreter. The `makeplots_matchpaper` script can be run as follows:

```
python makeplots_matchpaper.py
```

with outputs stored in the `plots` directory. Note that the `makeplots_matchpaper` is currently set up to produce one plot from histogram.dat files created using the `oxford_atlas_qcd` algorithm stored in `baseline_noPU_atlas_qcd`, and one plot from histogram.dat files created using the `oxford` algorithm stored in `baseline_noPU_oxford`. The two example 1D plots which are created are shown below:

`baseline_noPU_atlas_qcd_test.pdf`

![Example plot](plots/ExamplePlots/baseline_noPU_atlas_qcd_test.pdf)


`baseline_noPU_oxford_test.pdf` (note normalize has been set to True for this histogram, so the bin content is divided by the histogram integral)

![Example plot](plots/ExamplePlots/baseline_noPU_oxford_test.pdf)




