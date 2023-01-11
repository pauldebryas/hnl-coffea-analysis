# HNL analysis with coffea

Offline analysis of HNL signal, background and data.

##  Requirements

### Packages 

Here is a non-exhaustive list of the packages that I used:
- python3
- coffea
- akward array
- uproot
- numpy
- matplotlib 

For simplicity, you can use the same conda environment as I, which is saved in a .yml file.
For that:
- Make sure you have conda installed: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
- Install and activate the environment named "HNL" with the command:
```shell
conda env create -f HNL_environment.yml
conda activate HNL
```

### Inputs

Runs on pre-skimmed NanoAOD samples produced with https://github.com/cms-hnl/HNLTauPrompt
The list of NanoAOD samples that you want to use for the analysis must be specify in samples.py

##  Running the analysis

You can clone the repository with the folowing command:
```shell
git clone https://github.com/cms-hnl/hnl-coffea-tau
```

First you need to run `stitching/DY/stitching2D_DY.ipynb` and `stitching/WJets/stitching2D_WJets.ipynb` to compute and store stitching weights of Drell-Yan and W+Jets MC samples.

Then for the analysis, the files to be executed are `main/run_analysis.py` and `main/run_counters.py`. 
You can add --test option to make sure evrything works at first.
You need to specify tag that will be added to produced pkl files.

It saves the outcome in two pickle files located in the results/ directory:
- counter_{tag}.pkl which save the original sumw (before skimming) of MC samples (backgrounds and HNL signal) for scaling --> CountEvents.py file
- result_{tag}.pkl which save the histograms and sumw/nevent after each cut --> HNLAnalysis.py file

##  Results

###  main analysis
In the main/ folder
The files:
- plot_figures.ipynb allows you to see the histograms that are saved in the result_{tag}.pkl file and using counter_{tag}.pkl for scaling
- print_cutflow.ipynb allows you to analyse the cutflow also saved in the result_{tag}.pkl file and using counter_{tag}.pkl for scaling

###  DeepTau Comparaison 

In the same way, the DeepTauComparaison/ folder allows you to compare and see the improvement DeepTau2p5 brings to the analysis compare to 2p1

## Documentation
- Coffea: https://coffeateam.github.io/coffea/index.html
- pdgID of the particles: https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
- Akward array: https://awkward-array.readthedocs.io/en/latest/index.html
- Quick look at NanoAOD content: https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/data102X_doc.html#TrigObj\
- Correction for simulated tau object, as recommended by the TauPOG for UL2018; tau identification efficiency and fake rate scale factors (DeepTau2017v2p1VSe, DeepTau2017v2p1VSjet, DeepTau2017v2p1VSmu), tau energy scale, and tau trigger scale factors: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2\ (file with values https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/POG/TAU/2018_UL/tau.json.gz)
- Xsec: using DAS query for NanoAOD samples (https://cms-gen-dev.cern.ch/xsdb/) or using previous analysis (e.g https://github.com/hh-italian-group/hh-bbtautau/blob/master/Instruments/config/cross_section.cfg)
- Luminosity: In luminosity/run_Data there is the run number of all the runs in the data files (EGamma/SingleMuon/Tau) for 2018 and A/B/C/D area.
run2018_lumi.csv obtain at https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM give the corresponding luminosity of all those run number, so you can compute the exact luminosity for all data samples. 