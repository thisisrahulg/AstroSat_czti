#! /usr/bin/env python3
'''
Created by Rahul G ( on 22/04/2022 )

Run the code, provide new info, and copy and paste arranged_file.txt in czti archive page editor. back up of archive page is also created.

'''


from urllib.request import urlopen
from bs4 import BeautifulSoup
import pandas as pd
import re
import numpy as np
import sys
import Astrosat_Time as astime

url = 'http://astrosat.iucaa.in/czti/?q=grb'
html = urlopen(url)
soup = BeautifulSoup(html, 'lxml')



print('========= provide new information to be uploaded in GRB archive page ========== \n')

GRBname = str(input('GRB Name: '))
AstroSec= str(input('Astrosat seconds: '))
radec = str(input('RA, Dec :'))
thetaphi = str(input('Theta, phi : '))
tt90 = str(input('T90 : '))
GRB_linkname = str(input("GRB link's tag : "))

answer = input("Which lightcurves ? ( Enter 1 for czti, 2 for veto,3 for both) : ")
if answer == '1':
    tbin = str(input("  enter czt binning size [ 0.1 , 1 , 10 ] : "))
    cztlink = '<td class="rtecenter"><p><span style="font-size:12px;"><span style="font-family:times new roman,times,serif;"><a href="http://www.iucaa.in/~astrosat/czti_grb/{}/AS1CZT_{}_{}SECBIN_lightcurves.pdf" target="_blank">CZTI_lightcurve</a></span></span></p>'.format(GRB_linkname,GRB_linkname,tbin)
    finallink = cztlink

elif answer == '2':
    vtbin = str(input("  enter veto binning size [ 1 , 5 , 10 ]: "))
    vetolink = '<td class="rtecenter"><p><span style="font-size:12px;"><span style="font-family:times new roman,times,serif;"><a href="http://www.iucaa.in/~astrosat/czti_grb/{}/AS1CZT_{}_Veto_{}SECBIN_lightcurves.pdf" target="_blank">Veto_lightcurve</a></span></span></p>'.format(GRB_linkname,GRB_linkname,vtbin)
    finallink = vetolink

elif answer == '3':
    tbin = str(input("  enter czt binning size [ 0.1 , 1 , 10 ]: "))
    vtbin = str(input("  enter veto binning size [ 1 , 10 ]: "))
    cztlink = '<p><span style="font-size:12px;"><span style="font-family:times new roman,times,serif;"><a href="http://www.iucaa.in/~astrosat/czti_grb/{}/AS1CZT_{}_{}SECBIN_lightcurves.pdf" target="_blank">CZTI_lightcurve</a></span></span></p>'.format(GRB_linkname,GRB_linkname,tbin)
    vetolink = '<p><span style="font-size:12px;"><span style="font-family:times new roman,times,serif;"><a href="http://www.iucaa.in/~astrosat/czti_grb/{}/AS1CZT_{}_Veto_{}SECBIN_lightcurves.pdf" target="_blank">Veto_lightcurve</a></span></span></p>'.format(GRB_linkname,GRB_linkname,vtbin)
    finallink = '<td class="rtecenter">' + cztlink + vetolink
else:
    sys.exit("Give only mentioned values. Try again")

compton = input( 'Is compton curve is there? (y/n) : ')

if compton == 'y':
    complink = '<p><span style="font-size:12px;"><span style="font-family:times new roman,times,serif;"><a href="http://www.iucaa.in/~astrosat/czti_grb/{}/AS1CZT_{}_200_4_compton_events.pdf" target="_blank">Compton_lightcurve</a></span></span></p>'.format(GRB_linkname,GRB_linkname)
    finallink = finallink + complink + '</td>'

elif compton == 'n':
    finallink = finallink + '</td>'
    print('processing data ... \n')
else :
    sys.exit("Give only mentioned inputs. Try again")

print('\nadding new GRB info ...\n')

#----------------------- ADDING NEW DATA -----------------------------------------

initial = '<td class="rtecenter"><p><span style="font-size:12px;"><span style="font-family:times new roman,times,serif;">'
last = '</span></span></p></td>'

grb_name = initial + GRBname + last
astrosec = initial + AstroSec + last
rd = initial + radec + last
tp = initial + thetaphi + last
t90 = initial + tt90 + last
spectrum = finallink
rows=[]
tt =[]

new_row = [grb_name,astrosec,rd,tp,t90,spectrum]

print( '\nNew row created ... \n ')


table = soup.find_all('table')[0]

#-----------------CREATING BACKUP ---------------------------------------------------

print("\nCreating backup of archive data ... \n")

backup=open('archive_backup_utc.html', 'w')

backup.write(str(table))
backup.close()

print("\nBackup created ...\n")
#------------------------------------ARRANGING ARCHIVE PAGE---------------------------

for row in table.find_all('tr'):
    line = row.find_all('td')[1:]

    rows.append(line)
    time = line[1].text
    res= re.findall(r"[-+]?(?:\d*\.\d+|\d+)", time)
    tt.append(res)

time = AstroSec
result= re.findall(r"[-+]?(?:\d*\.\d+|\d+)", time)
tt.append(result)
tt.pop(0)
rows.append(new_row)
rows.pop(0)

finaltt = []


for i in range(len(tt)):
    finaltt.append(float(tt[i][0]))

sort_tt = sorted(finaltt)
sort_index = sorted(range(len(finaltt)), key=lambda k: finaltt[k])



arranged_table = []


for i,j in zip(sort_index,range(len(sort_index))):
    temp = rows[i]
    ttemp = []
    time_now=finaltt[i]
    for i in temp:
        ttemp.append(str(i))
    try:
        utc_time =  astime.convertAStoAll(time_now)[2]
        utc_time = utc_time[:-4]
    except:
        if time_now == 220934837.87 :
            utc_time = '2017-01-01 02:47:16'
        else:
            utc_time = 'NA'
    utc_string = '<td class="rtecenter" style="text-align: center;">\n<p><span style="font-size:12px;"><span style="font-family:times new roman,times,serif;">{}</span></span></p>\n</td>'.format(utc_time)
    ttemp.insert(2,utc_string)
    temp = " ".join(ttemp)
    temp =str(temp)

    final = "<tr><td class=\"rtecenter\"><p><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\">"+str(j+1)+"</span></span></p></td>"+ temp + "</tr>"
    arranged_table.append(final)


arranged_table.reverse()
joined = " ".join(arranged_table)



first = "<div class=\"rtecenter\">&nbsp;</div><div class=\"rtecenter\"><span style=\"color:#b22222;\"><span style=\"font-size:16px;\"><strong>&nbsp; &nbsp; &nbsp; &nbsp;&nbsp;&nbsp;ASTROSAT CZTI GRB ARCHIVE </strong></span></span></div><p>&nbsp;</p><table align=\"center\" border=\"1\" cellpadding=\"1\" cellspacing=\"1\" height=\"9519\" width=\"687\"><tbody><tr><td class=\"rtecenter\"><span style=\"font-size:16px;\"><strong>GRB No</strong></span></td><td class=\"rtecenter\"><span style=\"font-size:16px;\"><strong>GRB Name</strong></span></td><td class=\"rtecenter\"><p><span style=\"font-size:16px;\"><strong>Trigger time</strong></span></p><p><span style=\"font-size:16px;\"><strong>(Astrosat seconds)</strong></span></p></td><td class=\"rtecenter\"><p><span style=\"font-size:16px;\"><strong>RA, Dec</strong></span></p><p><span style=\"font-size:16px;\"><strong>(deg)</strong></span></p></td><td class=\"rtecenter\"><p><span style=\"font-size:16px;\"><strong>&nbsp;theta,&nbsp; phi&nbsp; </strong></span></p><p><span style=\"font-size:16px;\"><strong>(deg)</strong></span></p></td><td class=\"rtecenter\"><p><span style=\"font-size:16px;\"><strong>T90</strong></span></p><p><span style=\"font-size:16px;\"><strong>(sec)</strong></span></p></td><td class=\"rtecenter\"><span style=\"font-size:16px;\"><strong>Links</strong></span></td></tr>"

first = '<div class="rtecenter">&nbsp;</div><div class="rtecenter"><span style="color:#b22222;"><span style="font-size:16px;"><strong>&nbsp; &nbsp; &nbsp; &nbsp;&nbsp;&nbsp;ASTROSAT CZTI GRB ARCHIVE </strong></span></span></div><p>&nbsp;</p><table align="center" border="1" cellpadding="1" cellspacing="1" height="9519" width="687"><tbody><tr><td class="rtecenter"><span style="font-size:16px;"><strong>GRB No</strong></span></td><td class="rtecenter"><span style="font-size:16px;"><strong>GRB Name</strong></span></td><td class="rtecenter"><p><span style="font-size:16px;"><strong>Trigger time</strong></span></p><p><span style="font-size:16px;"><strong>(Astrosat seconds)</strong></span></p></td><td class="rtecenter"><p><span style="font-size:16px;"><strong>UTC Time</strong></span></p></td><td class="rtecenter"><p><span style="font-size:16px;"><strong>RA, Dec</strong></span></p><p><span style="font-size:16px;"><strong>(deg)</strong></span></p></td><td class="rtecenter"><p><span style="font-size:16px;"><strong>&nbsp;theta,&nbsp; phi&nbsp; </strong></span></p><p><span style="font-size:16px;"><strong>(deg)</strong></span></p></td><td class="rtecenter"><p><span style="font-size:16px;"><strong>T90</strong></span></p><p><span style="font-size:16px;"><strong>(sec)</strong></span></p></td><td class="rtecenter"><span style="font-size:16px;"><strong>Links</strong></span></td></tr>'

second ="</tbody></table><p><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\">&nbsp;</span></span></span></span></p><div class=\"rtecenter\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size: 16px;\"><span style=\"color: rgb(178, 34, 34);\"><strong>Description of Light Curve PDF</strong></span></span></span></span></span></span></div><div class=\"rtejustify\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\">&nbsp;</span></span></span></span></div><div class=\"rtejustify\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">The lightcurve PDF file for each GRB contains four pages for each of the four independent quadrants (A, B, C, D) of CZTI. The time of the GRB trigger is marked with a vertical dashed line. The plots for each quadrant are:</span></span></span></span><span></div><p><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\">&nbsp;</span></span></span></span></p><p class=\"rtejustify\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\"><strong>Page 1</strong>:&nbsp;</span><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">The upper left panel shows a &quot;spectrogram&quot; - counts as a function of energy and time. The time axis is in seconds since</span><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">&nbsp;</span><span class=\"Object\" id=\"OBJ_PREFIX_DWT52_com_zimbra_date\" role=\"link\" style=\"color: rgb(0, 90, 149); cursor: pointer; font-family: monospace; font-size: 14px;\">2010-01-01</span><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">&nbsp;</span><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">00:00:00 UT, while energy is in keV. By default, this is calculated with 20 keV bins for energy and 20 second bins for time. The lower left panel is the total lightcurve. The upper right panel is the total spectrum integrated over the entire observation. In other words, the upper left and lower right panels are projections of the upper left spectrogram. The lower right panel is a histogram of count rates. We can refer to the counts, spectrum and lightcurve as C(E,t), S(E) and L(t) respectively.</span></span></span></span></span></p><p><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\">&nbsp;</span></span></span></span></p><p class=\"rtejustify\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><strong style=\"font-size: 14px;\">Page 2:</strong><span style=\"font-size: 14px;\">&nbsp;</span><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">These plots are obtained by subtracting the mean CZTI spectrum from the spectrogram on page 1. In short, C_sub(E,t) = C(E,t) - S(E). This serves to highlight any transient activity. Panel details are as in page 1.</span></span></span></span></span></p><p class=\"rtejustify\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\">&nbsp;</span></span></span></span></p><p class=\"rtejustify\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><strong style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">Page 3:</strong><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">&nbsp;To estimate the significance of any transient activity, we take the mean-subtracted lightcurve at each energy, and divide it by the standard deviation at that energy. Mathematically, we have&nbsp;</span><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">C_norm(E,t) =&nbsp;</span><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">(C(E,t) - S(E)) / sigma(E)</span><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">. The other panels are projections of the 2d plot as on page 1.</span></span></span></span></span></p><p class=\"rtejustify\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\">&nbsp;</span></span></span></span></p><p class=\"rtejustify\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\"><strong>Page 4:&nbsp;</strong>This shows the total lightcurve (page 1 lower left), mean spectrum subtracted lightcurve (page 2 lower left), and the normalised lightcurve (page 3 lower left) on a single page.</span><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">&nbsp;</span><br style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\" />&nbsp;</span></span></span></span></p><p class=\"rtejustify\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><span style=\"font-size:12px;\"><span style=\"font-family:times new roman,times,serif;\"><strong><span style=\"color: rgb(51, 51, 51); font-family: monospace; font-size: 14px;\">Pages 5-8, 9-12 and 13-16 repeat these plots for quadrants B, C and D respectively.</span></strong></span></span></span></span></p>"

joined =  first + joined + second
d = open('arranged_file_utc.html','w')
d.write(joined)
d.close()

print('\nArchive page\'s updated data is created. Please copy and paste arranged_file.html to archive page editor.\n')


