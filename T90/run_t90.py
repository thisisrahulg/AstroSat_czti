import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
# parser.add_argument("infile", help="File with windows", type=str)
parser.add_argument("triggertime", help="Triggertime, name etc", type=str)
parser.add_argument("tbin", help="Binning to use", type=str)
parser.add_argument("interval", help="Number of bins for interval", type=int)
parser.add_argument("-q", help="Use only this quad", default=None, type=int)
parser.add_argument("-exclude", help="Exclude this quad", action='append', type=int)
parser.add_argument('-v','--veto',action='store_true',help="Specify this argument to use VETO lc")
# parser.add_argument("nruns", help="Number of runs", type=int)
args = parser.parse_args()


# times, tstarts, tstops = np.genfromtxt(args.infile, delimiter=',', unpack=True, skip_header=1)
# tbin = float(args.tbin)
time = args.triggertime
tbin = float(args.tbin)

quad_sel = ''
if args.q is not None:
    quads_to_not_use = [0,1,2,3]
    quads_to_not_use.remove(args.q)
    
    for i in quads_to_not_use:
        quad_sel += f' --exclude_quads {i}'

if args.exclude is not None:
    quad_sel = ''
    for i in args.exclude:
        quad_sel += f' --exclude_quads {i}'


def lc(time, tbin):
    orbit = os.listdir(f'{time}/czti/orbit/')[0]
    if(args.veto):
        return f'{time}/czti/orbit/{orbit}/{orbit}_veto_{tbin}'
    else:
        return f'{time}/czti/orbit/{orbit}/{orbit}_{tbin}_0'

def outfile(time,tbin):
	if(args.veto):
		return f'results/t90_{time}_{tbin}_veto'
	else:
		return f'results/t90_{time}_{tbin}_czt'

# def command_run_t90(time, tbin, tstart, tstop, nruns):
#     return f"python3 t90_sampling.py {lc(time, tbin)} {outfile(time, tbin)} --triggertime {time} --tstart {tstart} --tstop {tstop} --n_runs {nruns}"

def command_run_t90(time, tbin):
    return f"python3 t90_sampling.py {lc(time, tbin)} --outfile {outfile(time, tbin)} --triggertime {time} --interval {tbin*args.interval} {quad_sel}"

# for time, tstart, tstop in zip(times, tstarts, tstops):
#     print(f"Starting {time}")
#     command = command_run_t90(time, tbin, tstart, tstop, args.nruns)
#     print(command)
#     os.system(command)

command = command_run_t90(time, tbin)
os.system(command)
