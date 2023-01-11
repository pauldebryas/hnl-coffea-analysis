import sys
DIR_PATH = '/Users/debryas/Desktop/PhD_work/HNL_tau_analysis/hnl-coffea-analysis'
sys.path.append(DIR_PATH)

import numpy as np
import awkward as ak
from coffea import processor, hist

class Analysis_stitching2D_WJets(processor.ProcessorABC):
    def __init__(self, stitching_weights):
        ds_axis = hist.Cat("ds", "Primary dataset")
        acc_dict = {var: hist.Hist("Counts", ds_axis, axis) for var, axis in self.get_var_axis_pairs()}

        acc_dict[f'n_ev_all'] = processor.defaultdict_accumulator(int)
        acc_dict[f'sumw_all'] = processor.defaultdict_accumulator(float)

        self._accumulator = processor.dict_accumulator(acc_dict)

        self.stitching_weights = stitching_weights
    
    @property
    def accumulator(self):
        return self._accumulator

    @staticmethod
    def get_var_axis_pairs():
        leading_jet_pt_axis = hist.Bin("pt", "leading jet $p_{T}$ [GeV]", 15, 20, 150)
        NJet_axis = hist.Bin("NJet", "Number of jets", 12, 0, 12)
        Njets_axis = hist.Bin("Njets", "Number of LHE jets", 5, 0, 5)
        HT_axis = hist.Bin("HT", "LHE HT", 10, 0, 2600)

        v_a_pairs = [
            ('pt_leading_jet', leading_jet_pt_axis),
            ('NJet', NJet_axis),
            ('Njets', Njets_axis),
            ('HT', HT_axis),
        ]

        return v_a_pairs
    
    @staticmethod
    def get_NJets_bins():
        return [
            '0',
            '1',
            '2',
            '3',
            '4'
        ]

    @staticmethod
    def get_HT_bins():
        return [
            '0',
            '0to70',
            '70to100',
            '100to200',
            '200to400',
            '400to600',
            '600to800',
            '800to1200',
            '1200to2500',
            '2500toInf'
        ]

    # we will receive a NanoEvents instead of a coffea DataFrame
    def process(self, events):
        out = self.accumulator.identity()
        ds = events.metadata["dataset"] # dataset name
        #stitching
        WJets_samples = ['WJetsToLNu',
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

        if ds in WJets_samples:
            s = self.stitching_weights
            nHTbin = len(self.get_HT_bins())

            np.asarray(events.genWeight)[events.genWeight < 0] = -1.
            np.asarray(events.genWeight)[events.genWeight > 0] = 1.
            np.asarray(events.genWeight)[events.genWeight == 0] = 0.

            for NJets in self.get_NJets_bins():
                np.asarray(events.genWeight)[(events.LHE.HT == 0) & (events.LHE.Njets == int(NJets))] = events.genWeight[(events.LHE.HT == 0) & (events.LHE.Njets == int(NJets))]*s[nHTbin*int(NJets)+0]
                np.asarray(events.genWeight)[(events.LHE.HT > 0) & (events.LHE.HT < 70) & (events.LHE.Njets == int(NJets))] = events.genWeight[(events.LHE.HT > 0) & (events.LHE.HT < 70) & (events.LHE.Njets == int(NJets))]*s[nHTbin*int(NJets)+1]
                np.asarray(events.genWeight)[(events.LHE.HT >= 70) & (events.LHE.HT < 100) & (events.LHE.Njets == int(NJets))] = events.genWeight[(events.LHE.HT >= 70) & (events.LHE.HT < 100) & (events.LHE.Njets == int(NJets))]*s[nHTbin*int(NJets)+2]
                np.asarray(events.genWeight)[(events.LHE.HT >= 100) & (events.LHE.HT < 200) & (events.LHE.Njets == int(NJets))] = events.genWeight[(events.LHE.HT >= 100) & (events.LHE.HT < 200) & (events.LHE.Njets == int(NJets))]*s[nHTbin*int(NJets)+3]
                np.asarray(events.genWeight)[(events.LHE.HT >= 200) & (events.LHE.HT < 400) & (events.LHE.Njets == int(NJets))] = events.genWeight[(events.LHE.HT >= 200) & (events.LHE.HT < 400) & (events.LHE.Njets == int(NJets))]*s[nHTbin*int(NJets)+4]
                np.asarray(events.genWeight)[(events.LHE.HT >= 400) & (events.LHE.HT < 600) & (events.LHE.Njets == int(NJets))] = events.genWeight[(events.LHE.HT >= 400) & (events.LHE.HT < 600) & (events.LHE.Njets == int(NJets))]*s[nHTbin*int(NJets)+5]
                np.asarray(events.genWeight)[(events.LHE.HT >= 600) & (events.LHE.HT < 800) & (events.LHE.Njets == int(NJets))] = events.genWeight[(events.LHE.HT >= 600) & (events.LHE.HT < 800) & (events.LHE.Njets == int(NJets))]*s[nHTbin*int(NJets)+6]
                np.asarray(events.genWeight)[(events.LHE.HT >= 800) & (events.LHE.HT < 1200) & (events.LHE.Njets == int(NJets))] = events.genWeight[(events.LHE.HT >= 800) & (events.LHE.HT < 1200) & (events.LHE.Njets == int(NJets))]*s[nHTbin*int(NJets)+7]
                np.asarray(events.genWeight)[(events.LHE.HT >= 1200) & (events.LHE.HT < 2500) & (events.LHE.Njets == int(NJets))] = events.genWeight[(events.LHE.HT >= 1200) & (events.LHE.HT < 2500) & (events.LHE.Njets == int(NJets))]*s[nHTbin*int(NJets)+8]
                np.asarray(events.genWeight)[(events.LHE.HT >= 2500) & (events.LHE.Njets == int(NJets))] = events.genWeight[(events.LHE.HT >= 2500) & (events.LHE.Njets == int(NJets))]*s[nHTbin*int(NJets)+9]

        out['sumw_all'][ds] += ak.sum(events.genWeight)
        out['n_ev_all'][ds] += len(events)

        None_mask = ~ak.is_none( ak.flatten(events.GenJet.pt[ak.argmax(events.GenJet.pt, axis=-1, keepdims = True)]) )
        out[f'pt_leading_jet'].fill(ds=ds, pt=ak.flatten(events.GenJet.pt[ak.argmax(events.GenJet.pt, axis=-1, keepdims = True)])[None_mask], weight=events.genWeight[None_mask])
        out[f'NJet'].fill(ds=ds, NJet=ak.num(events.GenJet.pt), weight=events.genWeight)
        out[f'HT'].fill(ds=ds, HT=events.LHE.HT, weight=events.genWeight)
        out[f'Njets'].fill(ds=ds, Njets=events.LHE.Njets, weight=events.genWeight)

        return out

    def postprocess(self, accumulator):
        return accumulator