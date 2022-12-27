#! /usr/bin/env python3

import glob
import Astrosat_Time as astime
import sys

temp = glob.glob('*czt*')
cztfile = temp[0]

temp = glob.glob('*veto*')
vetofile = temp[0]

grbname =  input("Enter GRB name[GRB yyMMDD]: ")

title = '{}: AstroSat CZTI detection.\n\n'.format(grbname)
firstp = 'R. Gopalakrishnan (IUCAA), V. Prasad (IUCAA), G. Waratkar (IITB), A. Vibhute (IUCAA), V. Bhalerao (IITB), D. Bhattacharya (Ashoka University/IUCAA), A. R. Rao (IUCAA/TIFR), and S. Vadawale (PRL) report on behalf of the AstroSat CZTI collaboration:\n'
shortorlong = input('Is it short[1] or long[2]?: ')

if shortorlong == '1':
	secondp = 'Analysis of AstroSat CZTI data with the CIFT framework (Sharma et al., 2021, JApA, 42, 73) showed the detection of a short {} which was also detected by\n'.format(grbname)

elif shortorlong == '2':
	secondp = 'Analysis of AstroSat CZTI data with the CIFT framework (Sharma et al., 2021, JApA, 42, 73) showed the detection of a long {} which was also detected by\n'.format(grbname)

lastp = 'CZTI GRB detections are reported regularly on the payload site at http://astrosat.iucaa.in/czti/?q=grb. CZTI is built by a TIFR-led consortium of institutes across India, including VSSC, URSC, IUCAA, SAC, and PRL. The Indian Space Research Organisation funded, managed, and facilitated the project.\n'

lines = []
with open(cztfile) as file:
	for line in file:
		fields =  line.strip().split()
		lines.append(fields)


czttime = lines[5][-1]
czttime = astime.convertAStoAll(float(czttime))[2]



bkg = lines[7][4]
bkg1 = lines[7][5]
bkg2 = lines[7][6]
bkg1 = bkg1[2:]
bkg2 = bkg2[1:-1]
bkg = int(round(float(bkg),0))
bkg1 = int(round(float(bkg1),0))
bkg2 = int(round(float(bkg2),0))

total = lines[8][2]
total1 = lines[8][3]
total2 = lines[8][4]
total1 = total1[2:]
total2 = total2[1:-1]

total = int(round(float(total),0))
total1 = int(round(float(total1),0))
total2 = int(round(float(total2),0))

mean = lines[6][2]
mean1 = lines[6][3]
mean2 = lines[6][4]
mean1 = mean1[2:]
mean2 = mean2[1:-1]

mean = int(round(float(mean),0))
mean1 = int(round(float(mean1),0))
mean2 = int(round(float(mean2),0))

t90 = lines[9][1]
t901 = lines[9][2]
t902 = lines[9][3]
t901 = t901[2:]
t902 = t902[1:-1]

t90 = int(round(float(t90),0))
t901 = int(round(float(t901),0))
t902 = int(round(float(t902),0))

czt1 = 'The source was clearly detected in the 20-200 keV energy range.'

cztask = input('Is it multipeak in czt [y/n]: ')
quadccheck = input('Enter quadrants used in czti [2,3,4]: ')
if quadccheck == '2':
	quadc = 'two (out of four)'
elif quadccheck == '3':
	quadc = 'three (out of four)'
elif quadccheck == '4':
	quadc = 'all'
else:
	print('try again')
    sys.exit()

if cztask == 'y':
	czt2 = 'The light curve showed multiple peaks of emission with the strongest peak at {} UTC.'.format(czttime)
elif cztask == 'n':
	czt2 = ' The light curve peaks at {} UTC.'.format(czttime)

czt3 = ' The measured peak count rate associated with the burst is {} (+{}, -{}) counts/s above the background in the combined data of {} quadrants, with a total of {} (+{}, -{}) counts.'.format(bkg,bkg1,bkg2,quadc,total,total1,total2)
czt4 = ' The local mean background count rate was {} (+{}, -{}) counts/s. '.format(mean,mean1,mean2)
czt5 = 'Using cumulative rates, we measure a T90 of {} (+{}, -{}) s.'.format(t90,t901,t902)

veto = []
with open(vetofile) as v:
	for line in v:
		fields =  line.strip().split()
		veto.append(fields)

vetotime = veto[5][-1]
vetotime = astime.convertAStoAll(float(vetotime))[2]


bkg = veto[7][4]
bkg1 = veto[7][5]
bkg2 = veto[7][6]
bkg1 = bkg1[2:]
bkg2 = bkg2[1:-1]
bkg = int(round(float(bkg),0))
bkg1 = int(round(float(bkg1),0))
bkg2 = int(round(float(bkg2),0))

total = veto[8][2]
total1 = veto[8][3]
total2 = veto[8][4]
total1 = total1[2:]
total2 = total2[1:-1]

total = int(round(float(total),0))
total1 = int(round(float(total1),0))
total2 = int(round(float(total2),0))

mean = veto[6][2]
mean1 = veto[6][3]
mean2 = veto[6][4]
mean1 = mean1[2:]
mean2 = mean2[1:-1]

mean = int(round(float(mean),0))
mean1 = int(round(float(mean1),0))
mean2 = int(round(float(mean2),0))

t90 = veto[9][1]
t901 = veto[9][2]
t902 = veto[9][3]
t901 = t901[2:]
t902 = t902[1:-1]

t90 = int(round(float(t90),0))
t901 = int(round(float(t901),0))
t902 = int(round(float(t902),0))

veto1 = 'It was also clearly detected in the CsI anticoincidence (Veto) detector in the 100-500 keV energy range.'

vetocheck = input("Is veto multipeak or not? [y/n]: ")
quadvcheck = input('Enter quadrants used in veto  [2,3,4]: ')
if quadvcheck == '2':
	quadv = 'two (out of four)'
elif quadvcheck == '3':
	quadv = 'three (out of four)'
elif quadvcheck == '4':
	quadv = 'all'
else:
	print('try again')
    sys.exit()


if vetocheck == 'n':
	veto2 = ' The light curve peaks at {} UTC.'.format(vetotime)
else:
	veto2 = ' The light curve showed multiple peaks of emission with the strongest peak at {} UTC.'.format(vetotime)
veto3 = ' The measured peak count rate is {} (+{}, -{}) counts/s above the background in the combined Veto data of {} quadrants, with a total of {} (+{}, -{}) counts.'.format(bkg,bkg1,bkg2,quadv,total,total1,total2)

veto4 =  ' The local mean background count rate was {} (+{}, -{}) counts/s.'.format(mean,mean1,mean2)

veto5 = ' We measure a T90 of {} (+{}, -{}) s from the cumulative Veto light curve.'.format(t90,t901,t902)

with open('gcn_graft.txt', 'w') as g:
	g.write(title)
	g.write('\n\n')
	g.write(firstp)
	g.write('\n\n')
	g.write(secondp)
	g.write('\n\n')
	czt = czt1+czt2+czt3+czt4+czt5
	g.write(czt)
	g.write('\n\n')
	vet = veto1+veto2+veto3+veto4+veto5
	g.write(vet)
	g.write('\n\n')
	g.write(lastp)

print(" Check gcn_draft.txt for draft and add other detection details.")
