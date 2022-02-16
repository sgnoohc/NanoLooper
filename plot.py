#!/bin/env python

import plottery_wrapper as p

def plot(histname, xaxis_name):

    p.dump_plot(
        # reference frame (with lower value)
            fnames=["ttbaroutput.root"],
            # with larger value
            # data_fname="/home/users/joytzphysics/plots/VBSWWH_4p51.root",
            sig_fnames=["VBSWZHoutput.root"],
            legend_labels=["t#bar{t}"],
            signal_labels=["VBSWZH C_{2V}=3"],
            filter_pattern=histname,
            dogrep=False,
            extraoptions={
                "nbins": 60,
                "lumi_value": 137,
                #"ratio_range": [0., 100.],
                #"legend_datalabel": "VBSWWH C_{2V}=4p5",
                "legend_scalex":2.0,
                #"ratio_ndivisions":503,
                "xaxis_ndivisions":503,
                #"ratio_name": "C_{2V}=4p5 / C_{2V}=1",
                #"ratio_xaxis_title":xaxis_name,
                "print_yield":True,
                "yield_prec":4,
                },
    )

if __name__ == "__main__":

    hist_list = ("ptLep1",
                 "ptLep2",
                 "etaLep1",
                 "etaLep2",
                 "phiLep1",
                 "phiLep2",
                 "massDiLep",
                 "ptDiLep",
                 "dRLep",
                 "h_cutflow",
                 "ptVBFjet1",
                 "ptVBFjet2",
                 "etaVBFjet2",
                 "etaVBFjet1",
                 "WjetMass",
                 "VBFjetMass",
                 "ptWjet1",
                 "ptWjet2",
                 "dRleadingVBF",
                 "dRsubleadingVBF",
                 "etaVBFJet",
                 "ptHbb",
                 "etaHbb",
                 "massHbb",
                 "nFatjets",
                 "wjjScore",
                 "hbbScore",
                 "dRhiggs",
                 "ptWjj",
                 "etaWjj",
                 )
    for name in hist_list:
        plot(name, "GeV/c^2")

