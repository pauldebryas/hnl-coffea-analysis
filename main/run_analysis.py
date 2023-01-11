import sys
DIR_PATH = '/Users/debryas/Desktop/PhD_work/HNL_tau_analysis/hnl-coffea-analysis'
sys.path.append(DIR_PATH)

#!/usr/bin/env python
import pickle
import argparse
from coffea import processor
from coffea.nanoevents import NanoAODSchema
NanoAODSchema.warn_missing_crossrefs = True

from HNLAnalysis_ttm import HNLAnalysis_ttm
from HNLAnalysis_tte import HNLAnalysis_tte
from HNLAnalysis_tmm import HNLAnalysis_tmm
from HNLAnalysis_tee import HNLAnalysis_tee
from HNLAnalysis_tem_SS import HNLAnalysis_tem_SS
from HNLAnalysis_tem_OS import HNLAnalysis_tem_OS
#from HNLAnalysis_ttt import HNLAnalysis_ttt

from helpers import files_from_path
from samples.samples import signal_samples, Data_samples, MCbackground_samples

parser = argparse.ArgumentParser()
parser.add_argument('tag', help='tag that will be added to produced pkl files')
parser.add_argument('--test', action='store_true', default=False, help='test with a subset of the files')
parser.add_argument('-d', '--local_dir', default='/Users/debryas/cernbox/HNL/nanoV10/Run2_2018/')
parser.add_argument('-c', '--channel', default='ttm')
args = parser.parse_args()

local_dir = args.local_dir
test = args.test
tag = args.tag
channel = args.channel

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

if channel == 'tte':
    for region in regions:
        print('running region ' + region + ' for channel ' + channel)
        result = processor.run_uproot_job(
            samples,
            "Events",
            HNLAnalysis_tte(region),
            processor.iterative_executor, 
            {"schema": NanoAODSchema, 'workers': 6},
        )
        with open(f'{DIR_PATH}/results/result_{tag}_channel{channel}_region{region}.pkl', 'wb') as f:
            pickle.dump(result, f)

if channel == 'ttm':
    for region in regions:
        print('running region ' + region + ' for channel ' + channel)
        result = processor.run_uproot_job(
            samples,
            "Events",
            HNLAnalysis_ttm(region),
            processor.iterative_executor, 
            {"schema": NanoAODSchema, 'workers': 6},
        )
        with open(f'{DIR_PATH}/results/result_{tag}_channel{channel}_region{region}.pkl', 'wb') as f:
            pickle.dump(result, f)

if channel == 'tee':
    for region in regions:
        print('running region ' + region + ' for channel ' + channel)
        result = processor.run_uproot_job(
            samples,
            "Events",
            HNLAnalysis_tee(region),
            processor.iterative_executor, 
            {"schema": NanoAODSchema, 'workers': 6},
        )
        with open(f'{DIR_PATH}/results/result_{tag}_channel{channel}_region{region}.pkl', 'wb') as f:
            pickle.dump(result, f)

if channel == 'tmm':
    for region in regions:
        print('running region ' + region + ' for channel ' + channel)
        result = processor.run_uproot_job(
            samples,
            "Events",
            HNLAnalysis_tmm(region),
            processor.iterative_executor, 
            {"schema": NanoAODSchema, 'workers': 6},
        )
        with open(f'{DIR_PATH}/results/result_{tag}_channel{channel}_region{region}.pkl', 'wb') as f:
            pickle.dump(result, f)

if channel == 'tem_SS':
    for region in regions:
        print('running region ' + region + ' for channel ' + channel)
        result = processor.run_uproot_job(
            samples,
            "Events",
            HNLAnalysis_tem_SS(region),
            processor.iterative_executor,
            {"schema": NanoAODSchema, 'workers': 6},
        )
        with open(f'{DIR_PATH}/results/result_{tag}_channel{channel}_region{region}.pkl', 'wb') as f:
            pickle.dump(result, f)

if channel == 'tem_OS':
    for region in regions:
        print('running region ' + region + ' for channel ' + channel)
        result = processor.run_uproot_job(
            samples,
            "Events",
            HNLAnalysis_tem_OS(region),
            processor.iterative_executor,
            {"schema": NanoAODSchema, 'workers': 6},
        )
        with open(f'{DIR_PATH}/results/result_{tag}_channel{channel}_region{region}.pkl', 'wb') as f:
            pickle.dump(result, f)

if channel == 'ttt':
    for region in regions:
        print('running region ' + region + ' for channel ' + channel)
        result = processor.run_uproot_job(
            samples,
            "Events",
            HNLAnalysis_ttt(region),
            processor.iterative_executor,
            {"schema": NanoAODSchema, 'workers': 6},
        )
        with open(f'{DIR_PATH}/results/result_{tag}_channel{channel}_region{region}.pkl', 'wb') as f:
            pickle.dump(result, f)