#some useful functions for plots

import csv
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from coffea import hist
import warnings
import os

# for a given data samples, extract the good runs and the corresponding luminosity
def compute_lumi(ds):
    # for a given data sample, extract the good runs and the corresponding luminosity
    # good run for 2018 are stored in run2018_lumi.csv file
    with open('luminosity/run2018_lumi.csv', newline='') as csvfile:
        csv_reader = csv.reader(filter(lambda row: row[0]!='#', csvfile))
        run2018_goodrun = list(csv_reader)

    #store only information that we need: run number and luminosity
    run2018_run_lumi = []
    for i in range(len(run2018_goodrun)):
        run2018_run_lumi.append([run2018_goodrun[i][0][0:6],run2018_goodrun[i][5]])
    run2018_run_lumi = np.array(run2018_run_lumi).astype(float)

    #then found the run in the data file (stored in run_Data_ds.csv file)
    run_data = []
    with open('luminosity/run_Data/run_'+ds+'.csv', newline='') as csvfile:
        csv_reader = csv.reader(filter(lambda row: row[0]!='#', csvfile))
        run_data.append(list(csv_reader))
    run_data = np.concatenate(run_data).astype(float)

    # do the matching with the "good run" in run2018_lumi.csv
    run_lumi = []
    for i in range(len(run_data)):
        result = np.where(run2018_run_lumi[:,0] == run_data[i])
        if len(result[0]) == 1:
            index = result[0][0]
            run_lumi.append([run_data[i],run2018_run_lumi[index][1]])
        #if len(result[0]) == 0:
            #print("no good run found in run2018_lumi.csv for "+str(ds[:-1])+", run: "+str(run_data[i]))
        if len(result[0]) > 1:
            print("WARNING: "+str(ds[:-1])+", run: "+str(run_data[i])+" has multiple matching in run2018_lumi.csv")
    #array with good runs and their corresponding luminosity
    run_lumi = np.array(run_lumi).astype(float)
    #return the luminosity of ds
    return np.sum(run_lumi[:,1])*1000

def pt_ratio_plot(data_hist, BCK_histo, lumi, var, save_dir):

    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 18,
        'axes.labelsize': 18,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12
    })

    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(7,7),
        gridspec_kw={"height_ratios": (3, 1)},
        sharex=True
    )
    
    fig.subplots_adjust(hspace=.07)

    # Here is an example of setting up a color cycler to color the various fill patches
    # We get the colors from this useful utility: http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=6
    colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c']
    ax.set_prop_cycle(cycler(color=colors))

    fill_opts = {
        'edgecolor': (0,0,0,0.3),
        'alpha': 0.8
    }

    error_opts = {
        'label': 'Stat. Unc.',
        'hatch': '///',
        'facecolor': 'none',
        'edgecolor': (0,0,0,.5),
            'linewidth': 0
    }
    
    data_err_opts = {
        'linestyle': 'none',
        'marker': '.',
        'markersize': 10.,
        'color': 'k',
        'elinewidth': 1,
    }

    # now the data, setting clear=False to avoid overwriting the previous plot
    hist.plot1d(
        data_hist,
        overlay="ds",
        ax=ax,
        clear=False,
        #error_opts=data_err_opts
    )

    warnings.filterwarnings("ignore")

    hist.plot1d(
        #result[var][main_background],
        BCK_histo,
        overlay="ds",
        ax=ax,
        clear=False,
        stack=True,
        line_opts=None,
        fill_opts=fill_opts,
        #error_opts=error_opts
    )

    warnings.filterwarnings("ignore")

    ax.set_ylim(0, None)
    ax.set_xlim(0, 200)
    ax.set_xlabel(None)
    leg = ax.legend(prop={'size': 8})

    # now we build the ratio plot
    hist.plotratio(
        num=data_hist.sum("ds"),
        denom=BCK_histo.sum("ds"),
        ax=rax,
        error_opts=data_err_opts,
        denom_fill_opts={},
        guide_opts={},
        unc='num'
    )
    
    rax.set_ylabel('Ratio')
    rax.set_ylim(0,2)
    rax.set_xlim(0, 200)
            
    # add some labels
    #channel top left
    channel_topleft = plt.text(0., 1., var[3:],
                      fontsize=16,
                      horizontalalignment='left',
                      verticalalignment='bottom',
                      transform=ax.transAxes
                      )
    #lumi top right
    lumi_topright = plt.text(1., 1., str(round(lumi,2))+r" fb$^{-1}$ (13 TeV)",
                    fontsize=16,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    transform=ax.transAxes
                    )
    
    results_dir = os.path.join(save_dir, var[-3:]+'/') 
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    plt.savefig(results_dir+var+'_ratio.pdf',format='pdf')
    plt.show()