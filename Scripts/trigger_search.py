import sys,os
import argparse
parser = argparse.ArgumentParser(description='Script to plot lightcurve at a desired time')
parser.add_argument('file_name',help='file name of the pdf file',type=str)
parser.add_argument('observation',help='observation folder',type=str)
parser.add_argument('orbit',help='orbit number',type=str)
parser.add_argument('tmark',help='the mid time of lightcurve',type=float)
parser.add_argument('tbin',help='bin size [0.1,1,10], for veto [1,5,10]',type=float)
parser.add_argument('-v',"--veto",action='store_true',help='to plot veto lightcurve')

args = parser.parse_args()

path = '/data2/czti/level2/{obs}/czti/orbit/{orb}/modeM0/'.format(obs=args.observation,orb=args.orbit)
evtfile = ''
if os.path.exists(path):
        for fname in os.listdir(path):
                if fname.endswith('quad_clean.evt'):
                        evtfile = os.path.abspath(os.path.join(path,fname))
                        print('Event file using: {evtfile}'.format(evtfile=evtfile))
if evtfile == '':
        print('No event file found...\nexiting...')
        sys.exit()
if args.veto:
        if int(args.tbin) == 1:
                tmin = args.tmark - 150
                tmax = args.tmark + 150
        elif int(args.tbin) == 5:
                tmin = args.tmark - 300
                tmax = args.tmark + 300
        elif int(args.tbin) == 10:
                tmin = args.tmark - 500
                tmax = args.tmark + 500
else:
        if args.tbin == 0.1:
                tmin = args.tmark - 25
                tmax = args.tmark + 25
        elif int(args.tbin) == 1:
                tmin = args.tmark - 150
                tmax = args.tmark + 250
        elif int(args.tbin) == 10:
                tmin = args.tmark - 500
                tmax = args.tmark + 500

if args.veto:
        os.system('trans_det_V2.py {evtfile} {name}_{tbin}_veto --tmark {tmark} --tmin {tmin} --tmax {tmax} --tbin {tbin} --lc_type VETO'.format(name=args.file_name,evtfile=evtfile,tmark=args.tmark,tmin=tmin,tmax=tmax,tbin=args.tbin))
else:
        os.system('trans_det_V2.py {evtfile} {name}_{tbin}_czti --tmark {tmark} --tmin {tmin} --tmax {tmax} --tbin {tbin}'.format(name=args.file_name,evtfile=evtfile,tmark=args.tmark,tmin=tmin,tmax=tmax,tbin=args.tbin))




        
 




