import sys
DIR_PATH = '/Users/debryas/Desktop/PhD_work/HNL_tau_analysis/hnl-coffea-analysis'
sys.path.append(DIR_PATH)

import numpy as np
import awkward as ak
from coffea import processor, hist

class Analysis_stitching2D_DY_inclOnly(processor.ProcessorABC):
    def __init__(self):
        ds_axis = hist.Cat("ds", "Primary dataset")
        acc_dict = {var: hist.Hist("Counts", ds_axis, axis) for var, axis in self.get_var_axis_pairs()}

        acc_dict[f'n_ev_all'] = processor.defaultdict_accumulator(int)
        acc_dict[f'sumw_all'] = processor.defaultdict_accumulator(float)

        self._accumulator = processor.dict_accumulator(acc_dict)
    
    @property
    def accumulator(self):
        return self._accumulator

    @staticmethod
    def get_var_axis_pairs():
        leading_jet_pt_axis = hist.Bin("pt", "leading jet $p_{T}$ [GeV]", 15, 20, 150)
        NJet_axis = hist.Bin("NJet", "Number of jets", 12, 0, 12)
        NpNLO_axis = hist.Bin("NpNLO", "Number of LHE jets", 3, 0, 3)
        Vpt_axis = hist.Bin("Vpt", "LHE Z $p_{T}$ [GeV]", 30, 0, 800)

        v_a_pairs = [
            ('pt_leading_jet', leading_jet_pt_axis),
            ('NJet', NJet_axis),
            ('NpNLO', NpNLO_axis),
            ('Vpt', Vpt_axis),
        ]

        return v_a_pairs
    
    @staticmethod
    def get_NJets_bins():
        return [
            '0',
            '1',
            '2'
        ]

    @staticmethod
    def get_PtZ_bins():
        return [
            '0',
            '0to50',
            '50to100',
            '100to250',
            '250to400',
            '400to650',
            '650toInf'
        ]

    # we will receive a NanoEvents instead of a coffea DataFrame
    def process(self, events):
        out = self.accumulator.identity()
        ds = events.metadata["dataset"] # dataset name
        #stitching
        DY_samples =['DYJetsToLL_M-50',
                    'DYJetsToLL_0J',
                    'DYJetsToLL_1J',
                    'DYJetsToLL_2J',
                    'DYJetsToLL_LHEFilterPtZ-0To50',
                    'DYJetsToLL_LHEFilterPtZ-50To100',
                    'DYJetsToLL_LHEFilterPtZ-100To250',
                    'DYJetsToLL_LHEFilterPtZ-250To400',
                    'DYJetsToLL_LHEFilterPtZ-400To650',
                    'DYJetsToLL_LHEFilterPtZ-650ToInf']

        if ds in DY_samples:
            np.asarray(events.genWeight)[events.genWeight < 0] = -1.
            np.asarray(events.genWeight)[events.genWeight > 0] = 1.
            np.asarray(events.genWeight)[events.genWeight == 0] = 0.

        out['sumw_all'][ds] += ak.sum(events.genWeight)
        out['n_ev_all'][ds] += len(events)

        None_mask = ~ak.is_none( ak.flatten(events.GenJet.pt[ak.argmax(events.GenJet.pt, axis=-1, keepdims = True)]) )
        out[f'pt_leading_jet'].fill(ds=ds, pt=ak.flatten(events.GenJet.pt[ak.argmax(events.GenJet.pt, axis=-1, keepdims = True)])[None_mask], weight=events.genWeight[None_mask])
        out[f'NJet'].fill(ds=ds, NJet=ak.num(events.GenJet.pt), weight=events.genWeight)
        out[f'Vpt'].fill(ds=ds, Vpt=events.LHE.Vpt, weight=events.genWeight)
        out[f'NpNLO'].fill(ds=ds, NpNLO=events.LHE.NpNLO, weight=events.genWeight)
        
        return out

    def postprocess(self, accumulator):
        return accumulator