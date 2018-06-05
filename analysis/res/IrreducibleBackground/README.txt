The current samples being used are saved in a MATLAB Data file as: 
"sampleHighStats.mat"
And are not sorted by channel yet (done within the script)README:

The irreducible background estimation is done channel by channel in the following scripts:
BoostedEstimation
IntermediateEstimation
ResolvedEstimation

This involves creating mass templates for distributions in the leading and sub-leading invariant masses channel by channel. The analysis is done for the 4-tag samples. 

The likelihood fits are performed in:
BoostedLikelihood
IntermediateLikelihood
ResolvedLikelihood

And currently should be run after the appropriate estimation script. As new samples are used the start points would need to be changed, as would potentially the number of bins or fit function choices. 

Inputs are at `/eos/user/b/bstanisl/HH4b\ Pheno/Zach\ MPhys\ Input\ Files`

(The following scrips are just used as functions:
channelSort.m - Sort each sample by Resolved/Intermediate/Boosted
CrystalBall2.m - function of double tailed Crystal Ball for HH fits
HistoBins.m - Producing histograms in (mH0,mH1)
histogram1Dboth.m - Creates 1D histograms in mH0 and mH1
histogramPlots.m - Produces the plot of 2D histograms
histogramPlotsErr.m - Produces the plot of 2D histograms with: error/bin_Value.)
