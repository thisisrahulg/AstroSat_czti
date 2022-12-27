import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clip, sigma_clipped_stats
from scipy.optimize import curve_fit, OptimizeWarning
import warnings
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table

import argparse


# For plotting
ticksize = 16
titlesize = 20
labelsize = 18
plt.rcParams['xtick.labelsize'] = ticksize
plt.rcParams['ytick.labelsize'] = ticksize

def get_quad_data(filename, tmin, tmax):
    """
    Select trimmed light curve from .lc file
    """
    if "veto" in filename:
        tbl=Table.read(filename)
        this_r=np.array(tbl['orig_vetocts'])
        this_t=np.array(tbl['time'])
    else:
        data = fits.getdata(filename)
        this_t = data['time']
        this_r = data['rate']

    sel = (this_t<tmax) & (this_t>tmin)
    return this_t[sel], this_r[sel]

def get_total_data(lc, quads, tmin, tmax):
    """
    Get total light curve from quadrants in quads
    """
    t = []
    r = []
    for quad in quads:
        if "veto" in lc:
            filename=f'{lc}_Q{quad}_detrended.fits'
        else:
            filename = f'{lc}_Q{quad}.lc'
        this_t, this_r = get_quad_data(filename, tmin, tmax)
        t.append(this_t)
        r.append(this_r)

    return t[0], np.sum(r, axis=0)

def quadratic_trend(x, a, b, c):
    """
    Polynomial of Order 2 for fitting background.
    """
    return a*x**2 + b*x + c

def interp_cumsum(tbins, cumsum, quantile, multi=0):
    """
    Uses linear interpolation to find the tbin where cumsum reaches quantile.
    
    Parameters
    ----------
    tbins : array (n,)
        bin centers (for a light curve this is time)
    cumsum : array (n,)
        Cumulative Distribution. Nominally between 0 and 1. Noise can change this slightly
    quantile : float
        Value between 0 and 1. 
    multi : int
        If multiple crossings exist, choose the first if multi == 0 and last if multi==-1
        Code may break if different values given
    
    Returns
    -------
    t_quantile : float
        bin (time) where cumsum reaches quantile
    """
    index = np.where(np.diff(np.sign(cumsum - quantile)))[0][multi]
    t_quantile = np.interp(quantile, cumsum[index:index+2], tbins[index:index+2])
    return t_quantile

def get_trend(time, rate, grb_mask, return_all=False):
    """
    Get the background quadratic trend with outlier rejection in the background
    
    Parameters
    ----------
    time : array (n,)
    rate : array (n,)
    grb_mask : array (n,), dtype `bool`
        True for GRB window. Region excluded from background fit
    return_all : bool
        Changes the returned quantities (see Returns)
    
    Returns:
    -------
    if return_all:
        first_trend  : array (n,)
            The first fit for trend, without background outlier rejection
        ratio  : array (n,)
            Ratio of rate to trend. Used for outlier rejection
    final_trend : array (n,)
        Final trend with background outlier rejection
    bkg_mask : array (n,)
        Mask is True for background outliers
    """
    time = time-time[0]
    p0 = [0, 0, np.median(rate[~grb_mask])]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        popt, pcov = curve_fit(quadratic_trend, time[~grb_mask], rate[~grb_mask], p0)
    first_trend = quadratic_trend(time, *popt)
    ratio = rate/first_trend
    rat_clipped = sigma_clip(ratio[~grb_mask], sigma=3, maxiters=100)

    bkg_mask = np.zeros_like(grb_mask)
    bkg_mask[~grb_mask] = rat_clipped.mask
    p0 = [0, 0, np.median(rate[~(bkg_mask|grb_mask)])]

    with warnings.catch_warnings():
        warnings.filterwarnings("error")
        try:
            popt, pcov = curve_fit(quadratic_trend, time[~(bkg_mask|grb_mask)], rate[~(bkg_mask|grb_mask)], p0)
            final_trend = quadratic_trend(time, *popt)
        except:
            raise OptimizeWarning(" Non-Convergence")

    if return_all:
        return first_trend, ratio, final_trend, bkg_mask
    else:
        return final_trend, bkg_mask

def get_stats(time, rate, grb_mask, for_plot=False):
    """
    Get relevant statistics for the GRB.
    
    Parameters
    ----------
    time : array (n,)
    rate : array (n,)
    grb_mask : array (n,), dtype `bool`
        True for GRB window. Region excluded from background fit
    for_plot : bool
        Changes the returned quantities (see Returns)
    
    Returns
    -------
    if for_plot:
        cumsum : array (n,)
            Outlier Rejected Cumulative Light Curve (normalized)
        bkg_mask : array (n,) dtype bool
            True where bin is outlier. time[~bkg_mask] should be used for plotting vs cumsum
        t05 : float
            Time where cumsum reaches 0.05 (5%)
        t95 : float
            Time where cumsum reaches 0.95 (95%)
    else:
        bkg_rate : float
            Background Rate (counts/s)
        peak_rate : float
            Peak Rate (counts/s)
        tot_counts : float
            Total GRB Counts (counts)
        t90 : float
            T90 of GRB (s)
        cumsum : array (n,)
            Outlier Rejected Cumulative Light Curve
        bkg_mask : array (n,)
            Mask is True for background outliers
        t05 : float
            Time where cumsum reaches 0.05 (5%)
        t95 : float
            Time where cumsum reaches 0.95 (95%)
    Discussion
    ----------
    Varun Bhalerao:
        T05 = first time within the GRB window that you cross the 5% level
        T95 = last time within the GRB window that you cross the 95% line
    Gaurav Waratkar:
        Shouldn't this be the 'first time within the GRB window that you cross 95% line'?
        Since our GRB window selection can be subjective and if there is a grb window that contains a noise at the end, there can be places within this window where the 95% line is crossed multiple times. To be conservative, we were picking a window that is slightly bigger than the obvious GRB.
    Varun Bhalerao: 
        My rationale was that one will pick the GRB window to be "definitely the GRB". But your view seems more valid - we usually select a more conservative window such that I am 100% sure nothing outside this is the GRB.
        So LAST 5% crossing and FIRST 95% crossing are the defining points to be used.
    """
    try:
        final_trend, bkg_mask = get_trend(time, rate, grb_mask)
    except OptimizeWarning:
        return [np.nan]*4
    detrend = rate-final_trend

    bkg_rate = np.mean(final_trend[~(grb_mask|bkg_mask)])
    peak_rate = np.max(detrend[grb_mask])
    tot_counts = np.sum(detrend[grb_mask])*np.median(time[1:] - time[:-1])
    tpeak = time[grb_mask][np.argmax(detrend[grb_mask])]
    
    cumsum = np.cumsum(detrend[~bkg_mask])
    grb_start_ind = np.where(grb_mask[~bkg_mask])[0][0]
    grb_stop_ind = np.where(grb_mask[~bkg_mask])[0][-1]

    _, pre_level, _ = sigma_clipped_stats(cumsum[:grb_start_ind])
    normalized_cumsum = cumsum - pre_level
    _, post_level, _ = sigma_clipped_stats(normalized_cumsum[grb_stop_ind+1:])
    normalized_cumsum = normalized_cumsum/post_level

    try:
        t05 = interp_cumsum(time[grb_mask], normalized_cumsum[grb_mask[~bkg_mask]], 0.05, multi=-1)    
        t95 = interp_cumsum(time[grb_mask], normalized_cumsum[grb_mask[~bkg_mask]], 0.95, multi=0)
        t90 = t95-t05
    except IndexError:
        t05 = np.nan
        t95 = np.nan
        t90 = np.nan

    if for_plot:
        return normalized_cumsum, bkg_mask, t05, t95
    else:
        return bkg_rate, peak_rate, tot_counts, t90, tpeak, cumsum, bkg_mask, t05, t95

def sample(time, rate, grb_mask, n_runs=5000, seed=27071999):
    """
    Find GRB parameters for `n_runs` generated light-curve.
    
    Parameters
    ----------
    time : array (n,)
    rate : array (n,)
    grb_mask : array (n,), dtype `bool`
        True for GRB window. Region excluded from background fit
    n_runs : int, default : 5000
        Number of runs to use for generating distributions
    seed : int, default : 27071999
        Random seed (for replicability)
    Returns
    -------
    poisson_rate : array (n, n_runs)
        Poisson Generated Light Curves
    bkg_rate_samp : array (n_runs,)
        Sampled Background Rate (counts/s)
    peak_rate_samp : array (n_runs,)
        Sampled Peak Rate (counts/s)
    tot_counts_samp : array (n_runs,)
        Sampled Total GRB Counts (counts)
    t90_samp : array (n_runs,)
        Sampled T90 of GRB (s)
    cumsum : list (length n_runs)
        List of Sampled Outlier Rejected Cumulative Light Curve (ragged list)
    bkg_mask : array (n_runs, n)
        Array of Background Outlier Masks for each sampling. Mask is true for outlier
    t05_samp : array (n_runs,)
        Sampled T05 of GRB (s)
    t95_samp : array (n_runs,)
        Sampled T95 of GRB (s)
    """
    np.random.seed(seed)
    tbin = np.median(time[1:] - time[:-1])
    poisson_rate = np.random.poisson(lam=tbin*rate.reshape(-1,1), size=(rate.size, n_runs))/tbin
    
    bkg_rate_samp = []
    peak_rate_samp = []
    tot_counts_samp = []
    t90_samp = []
    cumsum_samp = []
    bkg_mask_samp = []
    t05_samp = []
    t95_samp = []
    for i in range(n_runs):
        this_bkg_rate, this_peak_rate, this_tot_counts, this_t90, this_tpeak, this_cumsum, this_bkg_mask, this_t05, this_t95 = get_stats(time, poisson_rate[:,i], grb_mask)
        bkg_rate_samp = np.append(bkg_rate_samp, this_bkg_rate)
        peak_rate_samp = np.append(peak_rate_samp, this_peak_rate)
        tot_counts_samp = np.append(tot_counts_samp, this_tot_counts)
        t90_samp = np.append(t90_samp, this_t90)
        cumsum_samp.append(this_cumsum)
        bkg_mask_samp.append(this_bkg_mask)
        t05_samp = np.append(t05_samp, this_t05)
        t95_samp = np.append(t95_samp, this_t95)
        
    return poisson_rate, bkg_rate_samp, peak_rate_samp, tot_counts_samp, t90_samp, cumsum_samp, np.array(bkg_mask_samp), t05_samp, t95_samp

def plot_lc(time, rate, grb_mask):
    "Plot the Light Curve with GRB Window and Outlier Points marked"
    first_trend, ratio, final_trend, bkg_mask = get_trend(time, rate, grb_mask, return_all=True)
    
    plt.figure(figsize=(10,6))
    plt.plot(time, rate, color='tab:blue', lw=2.5)
    plt.axvspan(*time[grb_mask][[0,-1]].T, 0, 1, alpha=0.2, color='tab:green')
    plt.plot(time, first_trend, color='tab:orange',label="First Trend")
    plt.plot(time, final_trend, color='tab:purple',label="Final Trend")
    plt.scatter(time[bkg_mask], rate[bkg_mask], color='tab:red', label="Outliers")
    plt.scatter(time[grb_mask], rate[grb_mask], color='tab:green', label="GRB")
    plt.xlabel("AstroSat Time (s)", fontsize=labelsize)
    plt.ylabel("Rate (s$^{-1}$)", fontsize=labelsize)
    plt.title("Light Curve with GRB Window", fontsize=titlesize)
    plt.xlim(time[0], time[-1])
    plt.grid(axis='y', zorder=-10)
    plt.legend(fontsize=ticksize)

def plot_lc_detrend(time, rate, grb_mask):
    first_trend, ratio, final_trend, bkg_mask = get_trend(time, rate, grb_mask, return_all=True)
    detrend = rate-final_trend
    plt.figure(figsize=(10,6))
    plt.plot(time, detrend, color='tab:blue', lw=2.5)
    plt.axvspan(*time[grb_mask][[0,-1]].T, 0, 1, alpha=0.2, color='tab:green')
    # plt.plot(time, first_trend, color='tab:orange',label="First Trend")
    # plt.plot(time, final_trend, color='tab:purple',label="Final Trend")
    plt.scatter(time[bkg_mask], detrend[bkg_mask], color='tab:red', label="Outliers")
    plt.scatter(time[grb_mask], detrend[grb_mask], color='tab:green', label="GRB")
    plt.xlabel("AstroSat Time (s)", fontsize=labelsize)
    plt.ylabel("Excess Rate (s$^{-1}$)", fontsize=labelsize)
    plt.title("Light Curve with GRB Window", fontsize=titlesize)
    plt.xlim(time[0], time[-1])
    plt.grid(axis='y', zorder=-10)
    plt.legend(fontsize=ticksize)
    

def plot_lc_cum(time, rate, grb_mask):
    "Plot the Cumulative LC with T05 and T95 masked"
    cumsum, bkg_mask, t05, t95 = get_stats(time, rate, grb_mask, for_plot=True)
    plt.figure(figsize=(10,6))
    plt.plot(time[~bkg_mask], cumsum, color='tab:blue', lw=2.5, zorder=0)

    [plt.vlines(txx, 0, 1, color='k', linestyles=':', lw=2, zorder=10) for txx in [t05, t95]]

    plt.hlines([0.05, 0.95], time[0], time[-1], color='k', linestyles='--')
    plt.hlines([0.,1], time[0], time[-1], color='k', linestyles='-')
    plt.axvspan(*time[grb_mask][[0,-1]].T, 0, 1, alpha=0.2, color='tab:green', zorder=-1)
    plt.xlim(time[0], time[-1])
    plt.xlabel("AstroSat Time (s)", fontsize=labelsize)
    plt.ylabel("Cumulative Rate (s$^{-1}$)", fontsize=labelsize)
    plt.title(f"Cumulative Light Curve for T90 calculation\nT90: {t95-t05:.2f} s", fontsize=titlesize)

def plot_lc_sample(time, poisson_rate, grb_mask):
    "Plot Poisson Sampled LCs"
    plt.figure(figsize=(10,6))
    plt.plot(time, np.median(poisson_rate, axis=1), color='tab:red', lw=1)
    plt.axvspan(*time[grb_mask][[0,-1]].T, 0, 1, alpha=0.2, color='tab:green')
    plt.fill_between(time, 
                     np.median(poisson_rate, axis=1)-3*np.std(poisson_rate, axis=1), 
                     np.median(poisson_rate, axis=1)+3*np.std(poisson_rate, axis=1), 
                     color='tab:red', alpha=0.5)
    plt.xlabel("AstroSat Time (s)", fontsize=labelsize)
    plt.ylabel("Rate (s$^{-1}$)", fontsize=labelsize)
    plt.xlim(time[0], time[-1])
    plt.suptitle("Poisson Sampled Light Curve", fontsize=titlesize)
    plt.title("Light regions are 3 sigma limits for each bin", fontsize=titlesize-2)
    plt.grid(axis='y', zorder=-10)

def plot_lc_cum_sample(time, cumsum, grb_mask, bkg_mask):
    "Plot Poisson Sampled Cumulative LCs"
    bin_cumsum = []
    for i in range(len(cumsum)):
        this_cumsum = np.zeros(time.size)
        this_cumsum[bkg_mask[i]] = np.nan
        this_cumsum[~bkg_mask[i]] = cumsum[i]
        bin_cumsum.append(this_cumsum)
    bin_cumsum = np.array(bin_cumsum)

    median = np.nanmedian(bin_cumsum, axis=0)
    std = np.nanstd(bin_cumsum, axis=0)

    plt.figure(figsize=(10,6))
    plt.plot(time, median, color='tab:red', lw=2)
    plt.axvspan(*time[grb_mask][[0,-1]].T, 0, 1, alpha=0.2, color='tab:green')
    plt.fill_between(time, 
                     median-3*std, 
                     median+3*std, 
                     color='tab:red', alpha=0.5)

    plt.hlines([0.05*median[-1], 0.95*median[-1]], time[0], time[-1], color='k', linestyles='--')
    plt.hlines([0.,median[-1]], time[0], time[-1], color='k', linestyles='-')

    plt.axvline(np.nanmedian(t05_samp), 0, 1, color='k', ls='--')
    plt.axvline(np.nanmedian(t05_samp) - 3* np.nanstd(t05_samp), 0, 1, color='k', ls=':', lw=1)
    plt.axvline(np.nanmedian(t05_samp) + 3* np.nanstd(t05_samp), 0, 1, color='k', ls=':', lw=1)

    plt.axvline(np.nanmedian(t95_samp), 0, 1, color='k', ls='--')
    plt.axvline(np.nanmedian(t95_samp) - 3* np.nanstd(t95_samp), 0, 1, color='k', ls=':', lw=1)
    plt.axvline(np.nanmedian(t95_samp) + 3* np.nanstd(t95_samp), 0, 1, color='k', ls=':', lw=1)

    plt.xlim(time[0], time[-1])
    plt.xlabel("AstroSat Time (s)", fontsize=labelsize)
    plt.ylabel("Cumulative Rate (s$^{-1}$)", fontsize=labelsize)
    plt.suptitle("Poisson Sampled Cumulative Light Curve", fontsize=titlesize)
    plt.title("Light regions are 3 sigma limits for each bin", fontsize=titlesize-2)
    plt.grid(axis='y', zorder=-10)

def plot_param_dist(bkg_rate_samp, peak_rate_samp, tot_counts_samp, t90_samp):
    "Plot the distributions of the parameters"
    plt.figure(figsize=(10,6))
    plt.subplot(221)
    plt.grid(axis='y', zorder=-1)
    plt.hist(bkg_rate_samp, bins=50, color='tab:red', alpha=0.75, zorder=2, label="Background Rate")
    plt.legend()

    plt.subplot(222)
    plt.grid(axis='y', zorder=-1)
    plt.hist(peak_rate_samp, bins=50, color='tab:red', alpha=0.75, zorder=2, label="Peak Rate")
    plt.legend()

    plt.subplot(223)
    plt.grid(axis='y', zorder=-1)
    plt.hist(tot_counts_samp, bins=50, color='tab:red', alpha=0.75, zorder=2, label="Total GRB Counts")
    plt.legend()

    plt.subplot(224)
    plt.grid(axis='y', zorder=-1)
    plt.hist(t90_samp[~np.isnan(t90_samp)], bins=50, color='tab:red', alpha=0.75, zorder=2, label="T90")
    plt.legend()
    plt.suptitle("Parameter Distributions", fontsize=titlesize)

def asymmetric_errors(original_value, array):
    "Statistics for each parameter. Histogram is generated twice. Can fix?"
    counts, bins = np.histogram(array[~np.isnan(array)], bins=50)
    bins = 0.5*(bins[1:] + bins[:-1])
    
    cumsum = np.cumsum(counts)/np.sum(counts)
    c05 = interp_cumsum(bins, cumsum, 0.05, multi=-1)
    c95 = interp_cumsum(bins, cumsum, 0.95, multi=0)
    return original_value-c05, c95-original_value

def onclick(event):
    
    time_ = event.xdata
    global start_stop
    start_stop.append(time_)
    
    # Disconnect after 2 clicks
    if len(start_stop) == 2:
        fig.canvas.mpl_disconnect(cid)
        plt.close()
        
    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("lc", help="Path to lc file, truncated at _Q*.lc", type=str)
    parser.add_argument("--outfile", help="Stem to prepend to output txt and pdf files", default=None, type=str)
    parser.add_argument("--tmin", help="Start time of light curve to process. Ignored if --triggertime provided", default=None, type=float)
    parser.add_argument("--tmax", help="End time of light curve to process. Ignored if --triggertime provided", default=None, type=float)
    parser.add_argument("--triggertime", help="Triggertime of GRB. Overrules --tmin, --tmax", default=None, type=float)
    parser.add_argument("--interval", help="Interval around triggertime to choose.", default=100, type=float)
    parser.add_argument("--tstart", help="Start time of GRB window", default=None, type=float)
    parser.add_argument("--tstop", help="End time of GRB window", default=None, type=float)
    parser.add_argument("--n_runs", help="Number of Poisson Generated Samples to use", default=5000, type=int)
    parser.add_argument("--exclude_quads", help="Quadrants to Exclude. 0:A, 1:B, 2:C, 3:D", action='append')
    parser.add_argument("--no_recenter", help="If provided, will not recenter to time of peak of GRB", dest="recenter", action="store_false")

    args = parser.parse_args()


    if args.triggertime is not None:
        tmin = args.triggertime - args.interval
        tmax = args.triggertime + args.interval
    elif args.tmin is not None and args.tmax is not None:
        tmin = args.tmin
        tmax = args.tmax
    else:
        raise NotImplementedError("Default times for tmin and tmax are not implemented. Please provide either --triggertime (and optionally --interval) or --tmin AND --tmax")
    
    if args.outfile is None:
        if args.triggertime is not None:
            args.outfile = f'{args.triggertime}'    
        else:
            args.outfile = f'{(tmin+tmax)/2}'

    quads = [0,1,2,3]
    if args.exclude_quads is not None:
        for i in args.exclude_quads:
            try:
                quads.remove(int(i))
            except ValueError:
                warnings.warn(f"{i} is not a valid quadrant. Skipping", RuntimeWarning)


    time, rate = get_total_data(args.lc, quads, tmin, tmax)

    if args.tstart is not None and args.tstop is not None:
        tstart = args.tstart
        tstop = args.tstop
    else:
        # raise NotImplementedError("Automated GRB window selection is not implemented. Please provide BOTH --tstart and --tstop")

        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.step(time, rate)
        plt.xlabel('Time (AstroSat seconds)')
        plt.ylabel('Count rate')
        plt.axvline([args.triggertime], color='k', linestyle='dashed', label='Marked times',alpha=0.7)
        ax.set_xlim(tmin, tmax) 
        start_stop=[]

        #Call onclick function
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        
        plt.show()
        
        tstart = start_stop[0]
        tstop = start_stop[1]






    grb_mask = (time<tstop)&(time>tstart) 

    original_bkg_rate, original_peak_rate, original_tot_counts, original_t90, original_tpeak, original_cumsum, original_bkg_mask, original_t05, original_t95 = get_stats(time, rate, grb_mask)

    if args.recenter:
        tmin = tmin - args.triggertime + original_tpeak
        tmax = tmax - args.triggertime + original_tpeak
        time, rate = get_total_data(args.lc, quads, tmin, tmax)
        grb_mask = (time<tstop)&(time>tstart) 
        original_bkg_rate, original_peak_rate, original_tot_counts, original_t90, original_tpeak, original_cumsum, original_bkg_mask, original_t05, original_t95 = get_stats(time, rate, grb_mask)

    poisson_rate, bkg_rate_samp, peak_rate_samp, tot_counts_samp, t90_samp, cumsum, bkg_mask, t05_samp, t95_samp = sample(time, rate, grb_mask, n_runs=args.n_runs)



    params = ['bkg_rate', 'peak_rate', 'tot_counts', 't90']
    param_name = ["Background Rate", "Peak Rate above Background", "Total Counts", "T90"]
    param_unit = ["counts/s", "counts/s", "counts", "s"]
    mean_fmts = ['.1f', '.1f', '.0f', '.2f']
    err_fmts = ['.1f', '.1f', '.0f', '.2f']

    
    out_txt = f"Quadrants Used: {quads}\n"
    out_txt += f"Time Interval Used : {tmin}\t{tmax}\n"
    out_txt += f"Number of Runs: {args.n_runs}\n"
    out_txt += f"Fraction of T90 measurements rejected = {np.sum(np.isnan(t90_samp))/t90_samp.size*100:.1f}% \n"
    out_txt += f"Time Binning: {np.median(time[1:]-time[:-1]):2.1e} s\nTime of Peak Rate: {original_tpeak}\n"
    for i, param in enumerate(params):
        original = eval(f'original_{param}')
        lower, upper = asymmetric_errors(eval(f"original_{param}"), eval(f"{param}_samp"))
        out_txt += f"{param_name[i]}:  {original:{mean_fmts[i]}} (+{upper:{err_fmts[i]}} -{lower:{err_fmts[i]}}) {param_unit[i]}\n"
        
    print(out_txt)
    
    text = str(round(original_t90,2))
    with open(f'{args.outfile}_{text}.txt', 'w') as f:
        f.write(out_txt)
        

    

    plotfile = PdfPages(f'{args.outfile}_{text}.pdf')

    plot_lc(time, rate, grb_mask)
    plotfile.savefig()
    plt.clf()

    plot_lc_detrend(time, rate, grb_mask)
    plotfile.savefig()
    plt.clf()

    plot_lc_cum(time, rate, grb_mask)
    plotfile.savefig()
    plt.clf()

    plot_lc_sample(time, poisson_rate, grb_mask)
    plotfile.savefig()
    plt.clf()

    plot_lc_cum_sample(time, cumsum, grb_mask, bkg_mask)
    plotfile.savefig()
    plt.clf()

    plot_param_dist(bkg_rate_samp, peak_rate_samp, tot_counts_samp, t90_samp)
    plotfile.savefig()
    plt.clf()

    plotfile.close()
