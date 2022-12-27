import numpy as np 
import matplotlib.pyplot as plt
from astropy.io import fits
import glob

def get_quad_data(filename, tmin, tmax):
    data = fits.getdata(filename)
    this_t = data['time']
    this_r = data['rate']
    sel = (this_t<tmax)&(this_t>tmin)
    return this_t[sel], this_r[sel]

def total_data(lc, tmin, tmax):
    t=[]
    r=[]
    for lc_name in lc:
        this_t, this_r = get_quad_data(lc_name, tmin, tmax)
        t.append(this_t)
        r.append(this_r)

    return t[0], np.sum(r, axis=0)

quads = [0,1,2,3]
superids = [262444262.03, 304699467.0, 319101992.0, 321149996.6, 322101631.0]

tbin=1.0

def onclick(event):
    
    time_ = event.xdata
    global start_stop
    
    start_stop.append(time_)
    
    # Disconnect after 2 clicks
    if len(start_stop) == 2:
        fig.canvas.mpl_disconnect(cid)
        plt.close(1)
        tstart = start_stop[0]
        tstop = start_stop[1]
        
        with open ('Lightcurves/output.txt', 'a') as f:
            f.write(",".join([str(triggertime), str(tstart), str(tstop)]) + '\n')
        
    return

for superid in superids[:]:
    lclist = glob.glob(f"Lightcurves/*{superid}/lc_*{tbin}_4_Q*.lc") 
    lclist.sort()

    triggertime = superid
    print(triggertime)
    tmin = triggertime-100
    tmax = triggertime+100
    lc=[]
    lc.extend([lclist[0], lclist[1], lclist[2], lclist[3]])
    print(lc)
    time, rate = total_data(lc, tmin, tmax)

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    plt.step(time, rate)
    plt.xlabel('Time (AstroSat seconds)')
    plt.ylabel('Count rate')
    plt.axvline([triggertime], color='k', linestyle='dashed', label='Marked times',alpha=0.7)
    ax.set_xlim(tmin, tmax) 
    start_stop=[]

    #Call onclick function
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    
    plt.show()
        
