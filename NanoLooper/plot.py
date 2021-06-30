#!/bin/env python

import plottery_wrapper as p

def plot(histname, xaxis_name):

    p.dump_plot(
	    # reference frame (with lower value)
            fnames=["/home/users/joytzphysics/NanoLooper/NanoLooper/wzh.root"],
    	    # with larger value
            # data_fname="/home/users/joytzphysics/plots/VBSWWH_4p51.root",
            # sig_fnames=["VBSWWH_11.root"],
            legend_labels=["VBSWZH"],
            filter_pattern=histname,
            dogrep=False,
            extraoptions={
                "nbins": 60,
                "lumi_value": 137,
                "ratio_range": [0., 100.],
                # "legend_datalabel": "VBSWWH C_{2V}=4p5",
                # "legend_scalex":2.0,
                "ratio_ndivisions":503,
                "xaxis_ndivisions":503,
                # "ratio_name": "C_{2V}=4p5 / C_{2V}=1",
                # "ratio_xaxis_title":xaxis_name,
                },
            )

if __name__ == "__main__":

	plot("massVBF", "mass [GeV/c^2]")
	plot("deltaEta", "[-10,10]")
	plot("ptW", "pt [GeV]")
	plot("ptH", "pt [GeV]")
	plot("ptZ", "pt [GeV]")

