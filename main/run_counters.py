import sys
DIR_PATH = '/Users/debryas/Desktop/PhD_work/HNL_tau_analysis/hnl-coffea-analysis'
sys.path.append(DIR_PATH)

#!/usr/bin/env python
import pickle
import argparse
from coffea import processor
from coffea.nanoevents import NanoAODSchema
NanoAODSchema.warn_missing_crossrefs = True

from CountEvents import CountEvents, CountEvents_stitched_samples

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

Drell_Yann_samples = [
    'DYJetsToLL_M-50',
    'DYJetsToLL_0J',
    'DYJetsToLL_1J',
    'DYJetsToLL_2J',
    'DYJetsToLL_LHEFilterPtZ-0To50',
    'DYJetsToLL_LHEFilterPtZ-50To100',
    'DYJetsToLL_LHEFilterPtZ-100To250',
    'DYJetsToLL_LHEFilterPtZ-250To400',
    'DYJetsToLL_LHEFilterPtZ-400To650',
    'DYJetsToLL_LHEFilterPtZ-650ToInf']

WJets_samples= [
    'WJetsToLNu',
    'W1JetsToLNu',
    'W2JetsToLNu',
    'W3JetsToLNu',
    'W4JetsToLNu',
    'WJetsToLNu_HT-70To100',
    'WJetsToLNu_HT-100To200',
    'WJetsToLNu_HT-200To400',
    'WJetsToLNu_HT-400To600',
    'WJetsToLNu_HT-600To800',
    'WJetsToLNu_HT-800To1200',
    'WJetsToLNu_HT-1200To2500',
    'WJetsToLNu_HT-2500ToInf']

Backgrounds_stitched = Drell_Yann_samples + WJets_samples

samples_stitched = {}
for element in Backgrounds_stitched:
    samples_stitched[element] = files_from_path(local_dir+element)

if test:
    for k, v in samples.items():
        samples[k] = v[:1]
    for k, v in samples_stitched.items():
        samples_stitched[k] = v[:1]

event_counter_NotSelected = processor.run_uproot_job(
    samples_stitched,
    'EventsNotSelected',
    CountEvents_stitched_samples(),
    processor.iterative_executor,
    {"schema": NanoAODSchema, 'workers': 6},
)

event_counter_Selected = processor.run_uproot_job(
    samples_stitched,
    'Events',
    CountEvents_stitched_samples(),
    processor.iterative_executor,
    {"schema": NanoAODSchema, 'workers': 6},
)

event_counter = processor.run_uproot_job(
    samples,
    'Runs',
    CountEvents(),
    processor.iterative_executor,
    {"schema": NanoAODSchema, 'workers': 6},
)

for sample in list(samples_stitched):
    event_counter['sumw'][sample] = event_counter_NotSelected['sumw_fixed'][sample] + event_counter_Selected['sumw_fixed'][sample]

with open(f'{DIR_PATH}/results/counter_{tag}.pkl', 'wb') as f:
    pickle.dump(event_counter, f)

