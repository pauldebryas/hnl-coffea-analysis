# Xsec from here: https://cms-gen-dev.cern.ch/xsdb/
xsecs = {
    #HNL
    'HNL100':0.0002046,
    'HNL125':0.0001166,
    'HNL150':6.446e-05,
    'HNL200':2.356e-05,
    'HNL250':1.073e-05,
    'HNL300':5.563e-06,
    'HNL350':3.177e-06,
    'HNL400':1.962e-06,
    'HNL450':1.259e-06,
    'HNL500':8.346e-07,
    'HNL600':4.045e-07,
    'HNL700':2.159e-07,
    'HNL800':1.209e-07,
    'HNL900':7.169e-08,
    'HNL1000':4.387e-08,
    #MCbackground
    #Drell-Yann
    'DY1Jets_To_LL_M-50':928.3,
    'DY2Jets_To_LL_M-50':293.6,
    'DY3Jets_To_LL_M-50':86.53,
    'DY4Jets_To_LL_M-50':41.28,
    'DYJets_To_LL_M-10to50':15890.0,
    'DYJets_To_LL_M-50-amcatnloFXFX':6077.22, #ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#DY_Z
    'DYJets_To_LL_M-50-madgraphMLM':5398.0,
    'DYJets_To_LL_M-50-madgraphMLM-ext':5321.0,
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
    'W2Jets_To_LNu':2809.0,
    'W3Jets_To_LNu':826.3,
    'W4Jets_To_LNu':544.3,
    'WJ1ets_To_LNu':8927.0,
    'WJets_To_LNu':61526.7, #ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#W_jets
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
    'WWW_4F':0.2086,
    'WWZ_4F_ext':0.1707,
    'WWZ_4F':0.1651,
    'WZZ_ext':0.05709,
    'WZZ':0.05565,
    'ZZZ_ext':0.01476,
    'ZZZ':0.01398,
    #Data
    #EGamma
    'Data_EGamma_2018':1.,
    #SingleMuon
    'Data_SingleMuon_2018':1.,
    #Tau
    'Data_Tau_2018':1.,
}
