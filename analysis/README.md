README
======

In this folder we have the main analysis code 

HH4b.cc

which takes a Les Houches parton level event file, showers it with
Pythia8 and then does the corresponding analysis, filling histograms,
creating plots etc.

It requires that you have installed:

* Pythia8 (I have tested up to v8.201):

* YODA 1.3.0 (https://yoda.hepforge.org)

* HepMC 2.06.08 or later (http://lcgapp.cern.ch/project/simu/HepMC/)

* Fastjet http://fastjet.fr/

Here one also needs to install fastjet contrib
http://fastjet.hepforge.org/contrib/
to access the Variable R jets

The config programs for each of these dependencies should be reachable
by the makefile. You should edit the Makefile samples path to the location of your samples.

************************************

### To compile and run the code
```
make
```

This builds the analysis code

```
./run.sh
```

This submits batch jobs to run the analysis code for all the analyses specified in ./src/HH4b.cc.
The code will output all results to `/res/<analysis>/<sample>`
where `<analysis>` specifies the name of the type of analysis, and `<sample>` refers to the source .lhe file.

Results are presented as YODA FLAT histograms for the histograms specified in the analysis.
Additionally the kinematics for input to the MVA are outputed to an "ntuple.dat" file, on
a sample by sample basis.

Note that by default the `OxfordAnalysis` is run (i.e. the code in `src/oxford.cc is used) and `OxfordAtlasQcdAnalysis` is run (i.e. the code in `src/oxford_atlas_qcd.cc` is used). This is specified in `src/run.cc` by the following lines:

```
analyses.push_back(new OxfordAnalysis(run, sample, subsample));

analyses.push_back(new OxfordAtlasQcdAnalysis(run, sample, subsample, 2));
analyses.push_back(new OxfordAtlasQcdAnalysis(run, sample, subsample, 3));
analyses.push_back(new OxfordAtlasQcdAnalysis(run, sample, subsample, 4));

```
Where the numbers 2, 3, 4 correspond to the number of b-tags. The `OxfordAnalysis` algorithm produces the outputs used to make previous paper plots, including the inputs to the mva code. The `OxfordAtlasQCDAnalysis` algorithm produces the outputs needed by the background estimation. Users can comment out either of these algorithms if the results from them are not needed for their particular study, but remember to `make` the code after making the changes.

************************************

### Writing analyses

An example barebones analysis class can be found in `/src/basic.cc` and `/inc/basic.h`.
This may be expanded upon to make a new analysis. For use in the HH4b main code,
the new analysis should be added as a module in the Makefile, and added to the
list of analyses in HH4b.cc

************************************

### Merging and plotting results

Read the `README.md` file in `res`.

************************************

