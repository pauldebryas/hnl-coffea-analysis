# Xsec from here: https://cms-gen-dev.cern.ch/xsdb/ [pb]
xsecs = {
    #HNL
    'HNL_tau_M-100':0.0002046,
    'HNL_tau_M-125':0.0001166,
    'HNL_tau_M-150':6.446e-05,
    'HNL_tau_M-200':2.356e-05,
    'HNL_tau_M-250':1.073e-05,
    'HNL_tau_M-300':5.563e-06,
    'HNL_tau_M-350':3.177e-06,
    'HNL_tau_M-400':1.962e-06,
    'HNL_tau_M-450':1.259e-06,
    'HNL_tau_M-500':8.346e-07,
    'HNL_tau_M-600':4.045e-07,
    'HNL_tau_M-700':2.159e-07,
    'HNL_tau_M-800':1.209e-07,
    'HNL_tau_M-900':7.169e-08,
    'HNL_tau_M-1000':4.387e-08,
    #MCbackground
    #Drell-Yann
    'DYJetsToLL_0J':5129.0,
    'DYJetsToLL_1J':951.5,
    'DYJetsToLL_2J':361.4,
    'DY1Jets_To_LL_M-50':928.3,
    'DY2Jets_To_LL_M-50':293.6,
    'DY3Jets_To_LL_M-50':86.53,
    'DY4Jets_To_LL_M-50':41.28,
    'DYJets_To_LL_M-10to50':15890.0,
    #'DYJets_To_LL_M-50-amcatnloFXFX':6077.22, #ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#DY_Z
    'DYJetsToLL_M-50': 6404.0,
    'DYJets_To_LL_M-50-madgraphMLM':5398.0,
    'DYJets_To_LL_M-50-madgraphMLM-ext':5321.0,
    'DYJetsToLL_LHEFilterPtZ-0To50':1485.0,
    'DYJetsToLL_LHEFilterPtZ-50To100':397.4,
    'DYJetsToLL_LHEFilterPtZ-100To250':97.2,
    'DYJetsToLL_LHEFilterPtZ-250To400':3.701,
    'DYJetsToLL_LHEFilterPtZ-400To650':0.5086,
    'DYJetsToLL_LHEFilterPtZ-650ToInf':0.04728,
    #Electroweak
    'EWK_WMinus2Jets_W_To_LNu_M-50':32.05,
    'EWK_WPlus2Jets_W_To_LNu_M-50':39.05,
    'EWK_Z2Jets_Z_To_LL_M-50':3.987, #ref: HTT AN-19-109
    # #SingleTop
    'ST_t-channel_antitop_4f_InclusiveDecays':69.09,
    'ST_t-channel_top_4f_InclusiveDecays':115.3,
    'ST_tW_antitop_5f_inclusiveDecays':34.97,
    'ST_tW_top_5f_inclusiveDecays':34.91,
    #TTbar
    'TT_To_2L2Nu':88.29, #ref: HTT AN-19-109
    'TT_To_Hadronic':687.1,
    'TT_To_SemiLeptonic':365.35, #ref: HTT AN-19-109
    #TT+bosons
    'TTWJets_To_LNu':0.2161,
    'TTWW':0.007003,
    'TTWZ':0.002453,
    'TTZ_To_LLNuNu_M-10':0.2439,
    'TTZZ':0.001386,
    # W+jets
    'W1JetsToLNu':8927.0,
    'W2JetsToLNu':2809.0,
    'W3JetsToLNu':826.3,
    'W4JetsToLNu':544.3,
    'WJetsToLNu':53870.0,
    #'WJetsToLNu':61526.7, #ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#W_jets
    'WJetsToLNu_HT-70To100':1264.0,
    'WJetsToLNu_HT-100To200':1256.0,
    'WJetsToLNu_HT-200To400':335.5,
    'WJetsToLNu_HT-400To600':45.25,
    'WJetsToLNu_HT-600To800':10.97,
    'WJetsToLNu_HT-800To1200':4.933,
    'WJetsToLNu_HT-1200To2500':1.16,
    'WJetsToLNu_HT-2500ToInf':0.008001,
     #DiBoson
    'WW':75.95,
    'WW_To_2L2Nu':11.09,
    'WZ':27.59,
    'WZ_To_2Q2L':6.419,
    'WZ_To_3LNu':5.213,
    'ZZ':12.17,
    'ZZ_To_2L2Nu':0.9738,
    'ZZ_To_2Q2L':2.214,
    'ZZ_To_4L':1.325,
    #TriBoson
    'WWW_4F_ext':0.2158,
    'WWW':0.2158,
    'WWZ_4F_ext':0.1707,
    'WWZ':0.1707,
    'WZZ_ext':0.05709,
    'WZZ':0.05565,
    'ZZZ_ext':0.01476,
    'ZZZ':0.01398,
    #QCD
    'QCD_Pt_15to30':1244000000.0,
    'QCD_Pt_30to50':106500000.0,
    'QCD_Pt_50to80':15700000.0,
    'QCD_Pt_80to120':2346000.0,
    'QCD_Pt_120to170':407700.0,
    'QCD_Pt_170to300':103700.0,
    'QCD_Pt_300to470':1244000000.0,
    #Data
    #EGamma
    'EGamma_2018':1.,
    #SingleMuon
    'SingleMuon_2018':1.,
    #Tau
    'Tau_2018':1.,
}
