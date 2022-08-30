#!/usr/bin/env python
import pickle
import argparse

from coffea import processor
from coffea.nanoevents import NanoAODSchema
NanoAODSchema.warn_missing_crossrefs = True

from HNLAnalysis import HNLAnalysis
from CountEvents import CountEvents

from helpers import files_from_path
from samples import signal_samples, Data_samples, MCbackground_samples

parser = argparse.ArgumentParser()
parser.add_argument('tag', help='tag that will be added to produced pkl files')
parser.add_argument('--test', action='store_true', default=False, help='test with a subset of the files')
parser.add_argument('-d', '--local_dir', default='/Users/debryas/cernbox/HNL/skimmed_samples/nanoAOD/2018/')
args = parser.parse_args()

local_dir = args.local_dir
test = args.test
tag = args.tag

samples = {}

for element in signal_samples:
    samples[element] = files_from_path(local_dir+'HNL_tau/'+signal_samples[element])

for element in MCbackground_samples:
    samples[element] = files_from_path(local_dir+'MC_background/'+MCbackground_samples[element])

for element in Data_samples:
    samples[element] = files_from_path(local_dir+'data/'+Data_samples[element])

if test:
    for k, v in samples.items():
        samples[k] = v[:1]

result = processor.run_uproot_job(
    samples,
    "Events",
    HNLAnalysis(),
    #processor.futures_executor,
    processor.iterative_executor, # may be better for debugging
    {"schema": NanoAODSchema, 'workers': 6},
)

event_counter = processor.run_uproot_job(
    samples,
    'Runs',
    CountEvents(),
    #processor.futures_executor,
    processor.iterative_executor, # may be better for debugging

    {"schema": NanoAODSchema, 'workers': 6},
)

with open(f'results/counter_{tag}.pkl', 'wb') as f:
    pickle.dump(event_counter, f)

with open(f'results/result_{tag}.pkl', 'wb') as f:
    pickle.dump(result, f)