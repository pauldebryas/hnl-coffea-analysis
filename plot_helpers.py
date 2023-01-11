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
    with open('/luminosity/run2018_lumi.csv', newline='') as csvfile:
        csv_reader = csv.reader(filter(lambda row: row[0]!='#', csvfile))
        run2018_goodrun = list(csv_reader)

    #store only information that we need: run number and luminosity
    run2018_run_lumi = []
    for i in range(len(run2018_goodrun)):
        run2018_run_lumi.append([run2018_goodrun[i][0][0:6],run2018_goodrun[i][5]])
    run2018_run_lumi = np.array(run2018_run_lumi).astype(float)

    #then found the run in the data file (stored in run_Data_ds.csv file)
    run_data = []
    with open('/luminosity/run_Data/run_'+ds+'.csv', newline='') as csvfile:
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
    run_lumi = np.array(run_lumi, dtype='float').astype(float)
    #return the luminosity of ds
    return np.sum(run_lumi[:,1])*1000 #convertion /fb --> /pb

def ratio_plot(data_hist, BCK_histo, xlim, lumi, var, save_dir):

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
    colors = ['#b2df8a','#1f78b4','#33a02c','#a6cee3','#e31a1c','#fb9a99']
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

    ymax_data = np.max(list(data_hist.values().values())[0])
    ax.set_ylim(0, ymax_data + 0.3*ymax_data)
    ax.set_xlim(0, xlim)
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
    rax.set_xlim(0, xlim)
            
    # add some labels
    #channel top left
    channel_topleft = plt.text(0., 1., 'CMS Preliminary',
                      fontsize=16,
                      horizontalalignment='left',
                      verticalalignment='bottom',
                      transform=ax.transAxes
                      )
    #lumi top right
    lumi_topright = plt.text(1., 1., str(round(lumi*1e-3,1))+r" fb$^{-1}$ (13 TeV)",
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

def ratio_plot_plus_signal(data_hist, BCK_histo, signal_histo, var, save_dir, xlim = None, ylim = None, luminosity = 59832.475347, hnl_list = ['HNL_tau_M-300']):

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
    colors = ['#8dd3c7','#984ea3','#b2df8a', '#fdb462','#a6cee3','#e31a1c','#fb9a99','#6a3d9a']
    ax.set_prop_cycle(cycler(color=colors))

    fill_opts = {
        'edgecolor': (0,0,0,0.3),
        'alpha': 0.8
    }

    error_option = {
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

    ymax_bck = np.max(list(BCK_histo.values().values())[0])
    hist.plot1d(
        BCK_histo,
        overlay="ds",
        ax=ax,
        clear=False,
        stack=True,
        line_opts=None,
        fill_opts=fill_opts,
        error_opts=error_option
    )

    warnings.filterwarnings("ignore")

    # now the data, setting clear=False to avoid overwriting the previous plot
    ymax_data = np.max(list(data_hist.values().values())[0])
    hist.plot1d(
        data_hist,
        overlay="ds",
        ax=ax,
        clear=False,
        error_opts=data_err_opts
    )

    #plot signal rescale
    ymax_signal = np.max(list(signal_histo.values().values())[0])
    scale_factor = (ymax_data / ymax_signal) * 0.9
    signal_histo.scale({s:scale_factor for s in hnl_list}, axis='ds')
    hist.plot1d(
        signal_histo,
        overlay="ds",
        ax=ax,
        clear=False,
        #error_opts=data_err_opts
    )

    warnings.filterwarnings("ignore")

    ax.set_ylim(0, ymax_data + 1.*ymax_data)
    if ylim != None:
        ax.set_ylim(0, ymax_data + ylim*ymax_data)

    if xlim != None:
        ax.set_xlim(xlim[0], xlim[1])
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
    
    rax.set_ylabel('Data/Bkg')
    rax.set_ylim(0,2)
    if xlim != None:
        rax.set_xlim(xlim[0], xlim[1])
            
    # add some labels
    #channel top left
    channel_topleft = plt.text(0., 1., 'CMS Preliminary',
                      fontsize=16,
                      horizontalalignment='left',
                      verticalalignment='bottom',
                      transform=ax.transAxes
                      )
    #lumi top right
    lumi_topright = plt.text(1., 1., str(round(luminosity*1e-3,1))+r" fb$^{-1}$ (13 TeV)",
                    fontsize=16,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    transform=ax.transAxes
                    )
    
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)

    plt.savefig(save_dir+var+'_ratio.pdf',format='pdf')
    plt.show()

def ratio_plot_DeepTau(data_hist, BCK_histo, var, save_dir, xlim = None, ylim = None, luminosity = 59832.475347):

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
    colors = ['#8dd3c7','#984ea3','#b2df8a', '#fdb462','#a6cee3','#e31a1c','#fb9a99','#6a3d9a']
    ax.set_prop_cycle(cycler(color=colors))

    fill_opts = {
        'edgecolor': (0,0,0,0.3),
        'alpha': 0.8
    }

    error_option = {
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

    ymax_bck = np.max(list(BCK_histo.values().values())[0])
    hist.plot1d(
        BCK_histo,
        overlay="ds",
        ax=ax,
        clear=False,
        stack=True,
        line_opts=None,
        fill_opts=fill_opts,
        error_opts=error_option
    )

    warnings.filterwarnings("ignore")

    # now the data, setting clear=False to avoid overwriting the previous plot
    ymax_data = np.max(list(data_hist.values().values())[0])
    hist.plot1d(
        data_hist,
        overlay="ds",
        ax=ax,
        clear=False,
        error_opts=data_err_opts
    )

    warnings.filterwarnings("ignore")

    ax.set_ylim(0, ymax_data + 1.*ymax_data)
    if ylim != None:
        ax.set_ylim(0, ymax_data + ylim*ymax_data)

    if xlim != None:
        ax.set_xlim(xlim[0], xlim[1])
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
    
    rax.set_ylabel('Data/Bkg')
    rax.set_ylim(0,2)
    if xlim != None:
        rax.set_xlim(xlim[0], xlim[1])
            
    # add some labels
    #channel top left
    channel_topleft = plt.text(0., 1., 'CMS Preliminary',
                      fontsize=16,
                      horizontalalignment='left',
                      verticalalignment='bottom',
                      transform=ax.transAxes
                      )
    #lumi top right
    lumi_topright = plt.text(1., 1., str(round(luminosity*1e-3,1))+r" fb$^{-1}$ (13 TeV)",
                    fontsize=16,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    transform=ax.transAxes
                    )
    
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)

    plt.savefig(save_dir+var+'_ratio.pdf',format='pdf')
    plt.show()

def ratio_plot_stitching(incl_hist, stitch_histo, var, save_dir, xlim = None, ylim = None, luminosity = 59832.475347):

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
    colors = ['#8dd3c7','#984ea3','#b2df8a', '#fdb462','#a6cee3','#e31a1c','#fb9a99','#6a3d9a']
    ax.set_prop_cycle(cycler(color=colors))

    fill_opts = {
        'edgecolor': (0,0,0,0.3),
        'alpha': 0.8
    }

    error_option = {
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

    ymax_bck = np.max(list(stitch_histo.values().values())[0])
    hist.plot1d(
        stitch_histo,
        overlay="ds",
        ax=ax,
        clear=False,
        stack=True,
        line_opts=None,
        fill_opts=fill_opts,
        error_opts=error_option
    )

    warnings.filterwarnings("ignore")

    # now the data, setting clear=False to avoid overwriting the previous plot
    ymax_data = np.max(list(incl_hist.values().values())[0])
    hist.plot1d(
        incl_hist,
        overlay="ds",
        ax=ax,
        clear=False,
        error_opts=data_err_opts
    )

    warnings.filterwarnings("ignore")

    ax.set_ylim(0, ymax_data + 1.*ymax_data)
    if ylim != None:
        ax.set_ylim(0, ymax_data + ylim*ymax_data)

    if xlim != None:
        ax.set_xlim(xlim[0], xlim[1])
    ax.set_xlabel(None)
    leg = ax.legend(prop={'size': 8})

    # now we build the ratio plot
    hist.plotratio(
        num=incl_hist.sum("ds"),
        denom=stitch_histo.sum("ds"),
        ax=rax,
        error_opts=data_err_opts,
        denom_fill_opts={},
        guide_opts={},
        unc='num'
    )
    
    rax.set_ylabel('incl/stitch')
    rax.set_ylim(0,2)
    if xlim != None:
        rax.set_xlim(xlim[0], xlim[1])
            
    # add some labels
    #channel top left
    channel_topleft = plt.text(0., 1., 'CMS Preliminary',
                      fontsize=16,
                      horizontalalignment='left',
                      verticalalignment='bottom',
                      transform=ax.transAxes
                      )
    #lumi top right
    lumi_topright = plt.text(1., 1., str(round(luminosity*1e-3,1))+r" fb$^{-1}$ (13 TeV)",
                    fontsize=16,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    transform=ax.transAxes
                    )
    
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)

    plt.savefig(save_dir+var+'_ratio_stitch_incl.pdf',format='pdf')
    plt.show()