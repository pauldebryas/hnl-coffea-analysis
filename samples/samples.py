#list of the orthogonal samples

signal_samples = {
    #'HNL_tau_M-85':'HNL_tau_M-85',
    'HNL_tau_M-100':'HNL_tau_M-100',
    #'HNL125':'HNL_tau_M-125',
    #'HNL150':'HNL_tau_M-150',
    'HNL_tau_M-200':'HNL_tau_M-200',
    #'HNL250':'HNL_tau_M-250',
    'HNL_tau_M-300':'HNL_tau_M-300',
    #'HNL350':'HNL_tau_M-350',
    'HNL_tau_M-400':'HNL_tau_M-400',
    #'HNL450':'HNL_tau_M-450',
    'HNL_tau_M-500':'HNL_tau_M-500',
    'HNL_tau_M-600':'HNL_tau_M-600',
    'HNL_tau_M-700':'HNL_tau_M-700',
    'HNL_tau_M-800':'HNL_tau_M-800',
    'HNL_tau_M-900':'HNL_tau_M-900',
    'HNL_tau_M-1000':'HNL_tau_M-1000',
}

MCbackground_samples = {
    #Drell-Yann
    'DYJetsToLL_M-50':'DYJetsToLL_M-50',
    'DYJetsToLL_0J':'DYJetsToLL_0J',
    'DYJetsToLL_1J':'DYJetsToLL_1J',
    'DYJetsToLL_2J':'DYJetsToLL_2J',
    'DYJetsToLL_LHEFilterPtZ-0To50':'DYJetsToLL_LHEFilterPtZ-0To50',
    'DYJetsToLL_LHEFilterPtZ-50To100':'DYJetsToLL_LHEFilterPtZ-50To100',
    'DYJetsToLL_LHEFilterPtZ-100To250':'DYJetsToLL_LHEFilterPtZ-100To250',
    'DYJetsToLL_LHEFilterPtZ-250To400':'DYJetsToLL_LHEFilterPtZ-250To400',
    'DYJetsToLL_LHEFilterPtZ-400To650':'DYJetsToLL_LHEFilterPtZ-400To650',
    'DYJetsToLL_LHEFilterPtZ-650ToInf':'DYJetsToLL_LHEFilterPtZ-650ToInf',
    #ELECRTOWEAK
    'EWK_WMinus2Jets_W_To_LNu_M-50':'EWKWMinus2Jets_WToLNu_M-50', 
    'EWK_WPlus2Jets_W_To_LNu_M-50':'EWKWPlus2Jets_WToLNu_M-50',
    'EWK_Z2Jets_Z_To_LL_M-50':'EWKZ2Jets_ZToLL_M-50',
    #singletop
    'ST_t-channel_antitop_4f_InclusiveDecays':'ST_t-channel_antitop_4f_InclusiveDecays',
    'ST_t-channel_top_4f_InclusiveDecays':'ST_t-channel_top_4f_InclusiveDecays',
    'ST_tW_antitop_5f_inclusiveDecays':'ST_tW_antitop_5f_inclusiveDecays',
    'ST_tW_top_5f_inclusiveDecays':'ST_tW_top_5f_inclusiveDecays',
    #TTbar
    'TT_To_2L2Nu':'TTTo2L2Nu',
    'TT_To_Hadronic':'TTToHadronic',
    'TT_To_SemiLeptonic':'TTToSemiLeptonic',
    #TT + bosons
    'TTWJets_To_LNu':'TTWJetsToLNu',
    'TTWW':'TTWW',
    'TTWZ':'TTWZ',
    'TTZ_To_LLNuNu_M-10':'TTZToLLNuNu_M-10',
    'TTZZ':'TTZZ',
    #W+jets
    'WJetsToLNu':'WJetsToLNu',
    'W1JetsToLNu':'W1JetsToLNu',
    'W2JetsToLNu':'W2JetsToLNu',
    'W3JetsToLNu':'W3JetsToLNu',
    'W4JetsToLNu':'W4JetsToLNu',
    'WJetsToLNu_HT-70To100':'WJetsToLNu_HT-70To100',
    'WJetsToLNu_HT-100To200':'WJetsToLNu_HT-100To200',
    'WJetsToLNu_HT-200To400':'WJetsToLNu_HT-200To400',
    'WJetsToLNu_HT-400To600':'WJetsToLNu_HT-400To600',
    'WJetsToLNu_HT-600To800':'WJetsToLNu_HT-600To800',
    'WJetsToLNu_HT-800To1200':'WJetsToLNu_HT-800To1200',
    'WJetsToLNu_HT-1200To2500':'WJetsToLNu_HT-1200To2500',
    'WJetsToLNu_HT-2500ToInf':'WJetsToLNu_HT-2500ToInf',
    #DiBoson
    'WW':'WW',
    'WZ':'WZ',
    'ZZ':'ZZ',
    #Tribosons
    'WWW':'WWW',
    'WWZ':'WWZ',
    'WZZ':'WZZ',
    'ZZZ':'ZZZ',
    #QCD
    #'QCD_Pt_15to30':'QCD_Pt_15to30',
    #'QCD_Pt_30to50':'QCD_Pt_30to50',
    #'QCD_Pt_50to80':'QCD_Pt_50to80',
    #'QCD_Pt_80to120':'QCD_Pt_80to120',
    #'QCD_Pt_120to170':'QCD_Pt_120to170',
    #'QCD_Pt_170to300':'QCD_Pt_170to300',
    #'QCD_Pt_300to470':'QCD_Pt_300to470'
}

Data_samples = {
    #EGamma
    'EGamma_2018A': 'EGamma_2018A',
    'EGamma_2018B':'EGamma_2018B',
    'EGamma_2018C':'EGamma_2018C',
    'EGamma_2018D':'EGamma_2018D',
    #SingleMuon
    'SingleMuon_2018A':'SingleMuon_2018A',
    'SingleMuon_2018B':'SingleMuon_2018B',
    'SingleMuon_2018C':'SingleMuon_2018C',
    'SingleMuon_2018D':'SingleMuon_2018D',
    #Tau
    'Tau_2018A':'Tau_2018A',
    'Tau_2018B':'Tau_2018B',
    'Tau_2018C':'Tau_2018C',
    'Tau_2018D':'Tau_2018D',
}