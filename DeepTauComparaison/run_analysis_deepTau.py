import sys
DIR_PATH = '/Users/debryas/Desktop/PhD_work/HNL_tau_analysis/hnl-coffea-analysis'
sys.path.append(DIR_PATH)

#!/usr/bin/env python
import pickle
import argparse
from coffea import processor
from coffea.nanoevents import NanoAODSchema
NanoAODSchema.warn_missing_crossrefs = True

from DeepTauComparaison.HNLAnalysis_ttm_new import HNLAnalysis_ttm_new 
from DeepTauComparaison.HNLAnalysis_ttm_old import HNLAnalysis_ttm_old
from DeepTauComparaison.HNLAnalysis_ttm_old_sftau import HNLAnalysis_ttm_old_sftau
from helpers import files_from_path
from samples.samples import signal_samples, Data_samples, MCbackground_samples

parser = argparse.ArgumentParser()
parser.add_argument('tag', help='tag that will be added to produced pkl files')
parser.add_argument('--test', action='store_true', default=False, help='test with a subset of the files')
parser.add_argument('-d', '--local_dir', default='/Users/debryas/cernbox/HNL/nanoV10/Run2_2018/')
args = parser.parse_args()

local_dir = args.local_dir
test = args.test
tag = args.tag

samples = {}

for element in signal_samples:
    samples[element] = files_from_path(local_dir+signal_samples[element])

for element in MCbackground_samples:
    samples[element] = files_from_path(local_dir+MCbackground_samples[element])

for element in Data_samples:
    samples[element] = files_from_path(local_dir+Data_samples[element])

regions = ['A','B','C','D']

if test:
    for k, v in samples.items():
        samples[k] = v[:1]
    regions = ['D']

for region in regions:
    print('running region ' + region + ' for channel ttm (old deep tau)')
    result = processor.run_uproot_job(
        samples,
        "Events",
        HNLAnalysis_ttm_old(region),
        processor.iterative_executor, 
        {"schema": NanoAODSchema, 'workers': 6},
    )

    with open(f'{DIR_PATH}/results/result_{tag}_old_region{region}.pkl', 'wb') as f:
        pickle.dump(result, f)


for region in regions:
    print('running region ' + region + ' for channel ttm (new deep tau)')
    result = processor.run_uproot_job(
        samples,
        "Events",
        HNLAnalysis_ttm_new(region),
        processor.iterative_executor, 
        {"schema": NanoAODSchema, 'workers': 6},
    )

    with open(f'{DIR_PATH}/results/result_{tag}_new_region{region}.pkl', 'wb') as f:
        pickle.dump(result, f)

for region in regions:
    print('running region ' + region + ' for channel ttm (old deep tau + tau scale factor)')
    result = processor.run_uproot_job(
        samples,
        "Events",
        HNLAnalysis_ttm_old_sftau(region),
        processor.iterative_executor, 
        {"schema": NanoAODSchema, 'workers': 6},
    )

    with open(f'{DIR_PATH}/results/result_{tag}_old_sftau_region{region}.pkl', 'wb') as f:
        pickle.dump(result, f)