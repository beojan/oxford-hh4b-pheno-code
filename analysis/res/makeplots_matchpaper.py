#!/usr/bin/env python2
from plotting_matchpaper import *

plot1D({'Signal': 'baseline_noPU_oxford/diHiggs/histo_m_H0_res_C1d.dat',
        'Background': 'baseline_noPU_oxford/background/histo_m_H0_res_C1d.dat'},normalize=True,
       xlabel=r'Invariant Mass of Leading Higgs Candidate (GeV)', title='Resolved category, no PU', output='test.pdf')
