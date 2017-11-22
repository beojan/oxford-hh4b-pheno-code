#!/usr/bin/env python2
from plotting_matchpaper import *

plot1D({'2b 2j': 'baseline_noPU_atlas_qcd/SHERPA_QCD2b2j/histo_pt_H0_res_RCO_4tag.dat',
        '4j': 'baseline_noPU_atlas_qcd/SHERPA_QCD4j/histo_pt_H0_res_RCO_4tag.dat'},
       xlabel=r'$p_{T}(h_{0})$ / GeV', title='4-tag All Regions', output='baseline_noPU_atlas_qcd_test.pdf')

plot1D({'Signal': 'baseline_noPU_oxford/diHiggs/histo_m_H0_res_C1d.dat',
       'Background': 'baseline_noPU_oxford/background/histo_m_H0_res_C1d.dat'},normalize=True,
      xlabel=r'Invariant Mass of Leading Higgs Candidate (GeV)', title='Resolved category, no PU', output='baseline_noPU_oxford_test.pdf')
