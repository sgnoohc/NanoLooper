#!/bin/env python

import plottery_wrapper as p

def plot(histname, xaxis_name):

    p.dump_plot(
	    # reference frame (with lower value)
            fnames=["/home/users/joytzphysics/NanoLooper/output.root"],
    	    # with larger value
            # data_fname="/home/users/joytzphysics/plots/VBSWWH_4p51.root",
            # sig_fnames=["VBSWWH_11.root"],
            legend_labels=["VBSWZH"],
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
		"print_yield":True
                },
    )

if __name__ == "__main__":

	hist_list = ( "h_gen_ptZ",
        "h_gen_ptW",
        "h_gen_ptH",
        "h_gen_deltaEta",
        "h_gen_massVBF",
	"h_gen_massbQsystem",
	"h_gen_ptb0",
	"h_gen_ptb1",
        "ptLep1",
        "ptLep2",
        "etaLep1",
        "etaLep2",
	"re_deltaEtaVBF",
        "phiLep1",
	"phiLep2",
	"massDiLep",
	"ptDiLep",
	"dRLep",
	"h_cutflow",
	"ptVBFjets",
	"etaVBFjets",
	"BjetMass",
	"BjetMass_400cut",
	"WjetMass",
	"VBFjetMass",
	"ptWjet1",
	"ptWjet2",
	"ptB1",
	"ptB2",
	"dRleadingVBF",
	"dRsubleadingbq"
    )
	for name in hist_list:
		plot(name, "GeV/c^2")

