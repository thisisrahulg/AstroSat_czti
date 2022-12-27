#! /usr/bin/env python3

import glob
import Astrosat_Time as astime
import sys
import re

temp = glob.glob('*czt*')
cztfile = temp[0]

temp = glob.glob('*veto*')
vetofile = temp[0]

grbname =  input("Enter GRB name[GRB yyMMDD]: ")

title = '{}: AstroSat CZTI detection.\n\n'.format(grbname)
firstp = 'R. Gopalakrishnan (IUCAA), V. Prasad (IUCAA), G. Waratkar (IITB), A. Suresh (IITB), A. Vibhute (IUCAA), V. Bhalerao (IITB), D. Bhattacharya (Ashoka University/IUCAA), A. R. Rao (IUCAA/TIFR), and S. Vadawale (PRL) report on behalf of the AstroSat CZTI collaboration:\n'


lastp = 'CZTI GRB detections are reported regularly on the payload site at http://astrosat.iucaa.in/czti/?q=grb. CZTI is built by a TIFR-led consortium of institutes across India, including VSSC, URSC, IUCAA, SAC, and PRL. The Indian Space Research Organisation funded, managed, and facilitated the project.\n\nLinks:\n------\n[1] http://astrosat.iucaa.in/czti/?q=grb\n'

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

checkt90 = float(t90) 

t90 = int(round(float(t90),0))
t901 = int(round(float(t901),0))
t902 = int(round(float(t902),0))

czt_quadrants = []
for q in range(2, len(lines[0])):
	quadr = lines[0][q]
	res = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", quadr)[0]
	if res == '0':
		czt_quadrants.append('A')
	elif res == '1':
		czt_quadrants.append('B')
	elif res == '2':
		czt_quadrants.append('C')
	elif res == '3':
		czt_quadrants.append('D')

string = ', '.join(czt_quadrants)
c0 = 'Quadrants Used: {}'.format(string)
c1 = 'Time Binning: {} s'.format(float(lines[4][2]))
c2 = 'Time of Peak Rate: {}'.format(lines[5][-1])
c3 = 'Background Rate: {} (+{}, -{}) counts/s'.format(mean, mean1, mean2)
c4 = 'Peak Rate above Background: {} (+{}, -{}) counts/s'.format(bkg, bkg1, bkg2)
c5 = 'Total Counts: {} (+{}, -{}) counts'.format(total, total1, total2)
c6 = 'T90: {} (+{}, -{}) s'.format(t90, t901, t902)

param_czt = [c0, c1, c2, c3, c4, c5, c6]

if checkt90 < 2.0:
	secondp = 'Analysis of AstroSat CZTI data with the CIFT framework (Sharma et al., 2021, JApA, 42, 73) showed the detection of a short {} which was also detected by\n'.format(grbname)

else:
	secondp = 'Analysis of AstroSat CZTI data with the CIFT framework (Sharma et al., 2021, JApA, 42, 73) showed the detection of a long {} which was also detected by\n'.format(grbname)

czt1 = 'The source was clearly detected in the 20-200 keV energy range.'

cztask = input('Is it multipeak in czt? [y/n]: ')
if len(lines[0]) == 4:
	quadc = 'two (out of four)'
elif len(lines[0]) == 5:
	quadc = 'three (out of four)'
elif len(lines[0]) == 6:
	quadc = 'all'
else:
	print('try again')
	sys.exit()

if cztask == 'y':
	czt2 = ' The light curve showed multiple peaks of emission with the strongest peak at {} UTC.'.format(czttime)
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

veto_quadrants = []
for q in range(2, len(veto[0])):
	quadr = veto[0][q]
	res = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", quadr)[0]
	if res == '0':
		veto_quadrants.append('A')
	elif res == '1':
		veto_quadrants.append('B')
	elif res == '2':
		veto_quadrants.append('C')
	elif res == '3':
		veto_quadrants.append('D')

string = ', '.join(veto_quadrants)
v0 = 'Quadrants Used: {}'.format(string)
v1 = 'Time Binning: {} s'.format(float(veto[4][2])) 
v2 = 'Time of Peak Rate: {}'.format(veto[5][-1])
v3 = 'Background Rate: {} (+{}, -{}) counts/s'.format(mean, mean1, mean2)
v4 = 'Peak Rate above Background: {} (+{}, -{}) counts/s'.format(bkg, bkg1, bkg2)
v5 = 'Total Counts: {} (+{}, -{}) counts'.format(total, total1, total2)
v6 = 'T90: {} (+{}, -{}) s'.format(t90, t901, t902)

param_veto = [v0, v1, v2, v3, v4 , v5, v6]

veto1 = 'It was also clearly detected in the CsI anticoincidence (Veto) detector in the 100-500 keV energy range.'

vetocheck = input("Is it multipeak in veto? [y/n]: ")

if len(veto[0]) == 4:
	quadv = 'two (out of four)'
elif len(veto[0]) == 5:
	quadv = 'three (out of four)'
elif len(veto[0]) == 6:
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

with open('gcn_draft.txt', 'w') as g:
	g.write('Parameters of the GRB (from CZT):\n\n')
	for i in param_czt:
		g.write(i)
		g.write('\n')
	g.write('\n\n')
	g.write('Parameters of the GRB (from Veto):\n\n')
	for i in param_veto:
		g.write(i)
		g.write('\n')
	g.write('\n\n')
	g.write('GCN Draft:\n\n')
	g.write(title)
	g.write(firstp)
	g.write('\n\n')
	g.write(secondp)
	g.write('\n')
	czt = czt1+czt2+czt3+czt4+czt5
	g.write(czt)
	g.write('\n\n')
	vet = veto1+veto2+veto3+veto4+veto5
	g.write(vet)
	g.write('\n\n')
	g.write(lastp)
	g.write('\n\n')
	g.write('Thanks,\nRahul\n\n')
	
print("\n Check gcn_draft.txt for draft and add other detection details.\n\n>>>>>>>>>>>>>>>>>>>> ADD OTHER DETECTIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")
