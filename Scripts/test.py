#!/usr/bin/env python2.7

"""
    Astrosat_Time.py
    Aim: This is a time conversion routine for Astrosat CZTI Instrument, It can be used to convert Astrosat Seconds (since epoch, 2010.0 UTC) to ISO 8601 date,
    Calender Date, Year and DayNumber, Julian Day and Modified Julian Day formats and Vice Versa.
    
    Type Astrosat_Time.py -h to get help

"""

#------------------------------------------------------------------------


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from urllib.request import urlopen
from bs4 import BeautifulSoup
import re
import os




start = int(input('Enter starting row: '))
stop = int(input('Enter last row: '))
print('\n--------------------------CALCULATING NOISE RATIO AND TELEMETRIC ERROR--------------------------\n')

file = open('all_noise_and_tele_error.dat','w')
file.write("Hi,\n \n")
file.write('---------------------------------------------\n')
file.write('                  Noise Ratio   \n')
file.write('---------------------------------------------\n')
file.write('ObsID   OrbitID    Quadrant    DetID   PixID   \n')
file.write('---------------------------------------------\n') 


   
url = 'http://www.iucaa.in/~astrosat/czti_dqr/'
html = urlopen(url)
s = BeautifulSoup(html, 'lxml')


links = s.find_all('a',href = True, limit = stop+6)

def mode(value):
	cleanv=[]
	for i in range(2,len(value)-1,4):
		cleanv.append(value[i].text)
		cleanv.append(value[i+1].text)
	error = []	
	for i in range(1,len(cleanv),2):
		p = float(cleanv[i])
		pp = float(cleanv[i-1])
		err = p/pp
		err = err*100.0
		err = round(err,6)
		error.append(err)
	return (error)
	
			
		
	

for l in range(stop+5,start+4,-1):
	
	
	strip = links[l].get('href')
	link = 'https://www.iucaa.in/~astrosat/czti_dqr/' + strip
	html_ = urlopen(link)

	soup = BeautifulSoup(html_, 'lxml')
	
	info = soup.find_all('li')[4]

	obsid = info.text[-5:-1]
	

	orb = soup.find_all('title')[0]
	orbid = orb.text[-5:]
	print('Checking ObsID:', obsid ,	' OrbID:', orbid)
	
	table = soup.find_all('table')
    
	table1 = table[1]
	mode1 = table1.find_all('tr')[0]
	val1 = table1.find_all('td')
	error1 = mode(val1)
		
	table2 = table[2]
	mode2 = table2.find_all('tr')[0]
	val2 = table2.find_all('td')
	error2 = mode(val2)
    
	t1 = 4
	t2 = 5    
	flag = 0 
    
	for i in range(len(error1)):
		if float(error1[i]) >= 5.0 or float(error2[i]) >= 5.0 :
		 
			flag = 1  
	if len(table) == 10 :
            
		t1 = 5
		t2 = 6

		
		table3 = table[3]
		mode3 = table3.find_all('tr')[0]
		val3 = table3.find_all('td')
		error3 = mode(val3) 
		for i in range(len(error1)):
			if float(error3[i]) >=5.0 :
				flag = 1
    	
	


    
	if flag == 1 :
		print (' Telemetry error is present in {} ,  orbit {}'.format(info.text[:-1],orbid))

		current_dir = os.getcwd()
		final_dir = os.path.join(current_dir, r'tele_error')
		if not os.path.exists(final_dir):
			os.makedirs(final_dir)
		
		file_name = '/tele_error_{}_{}.txt'.format(obsid,orbid)
		file_name = final_dir + file_name
		with open(file_name,'w') as er :
			er.write(' Telemetry error of {}, OrbitID : {}\n\n'.format(info.text[:-1],orbid))
			

			er.write(mode1.text +'\n' )
			er.write('-------------------------\n')	
			letter = ['A','B','C','D']
			
			for i in range(4):
				temp = ' {}    {} \n'.format(letter[i],error1[i])
				er.write(temp)
			
			er.write('\n')
			er.write(mode2.text +'\n' )
			er.write('-------------------------\n')			
				
			for i in range(4):
				temp = ' {}    {} \n'.format(letter[i],error2[i])
				er.write(temp)		

			
			if len(table) == 10 :
				er.write('\n')
				er.write(mode3.text +'\n' )
				er.write('-------------------------\n')
				for i in range(4):
					temp = ' {}    {} \n'.format(letter[i],error3[i])
					er.write(temp)

	

    	
	noise = soup.find_all('table')[t1]
	pixel = soup.find_all('table')[t2]

	detid = pixel.find_all('td')
	pixid = pixel.find_all('td')

	detid_a = detid[0].text
	pixid_a = pixid[1].text

	detid_b = detid[21].text
	pixid_b = pixid[22].text

	detid_c = detid[42].text
	pixid_c = pixid[43].text

	detid_d = detid[63].text
	pixid_d = pixid[64].text


	n1 = noise.find_all('td',class_='important')[0].text
	n2 = noise.find_all('td',class_='important')[1].text
	n3 = noise.find_all('td',class_='important')[2].text
	n4 = noise.find_all('td',class_='important')[3].text

	n1_ = float(n1[:-1])
	n2_ = float(n2[:-1])
	n3_ = float(n3[:-1])
	n4_ = float(n4[:-1])
	a = 'A = ' + n1
	b = 'B = ' + n2
	c = 'C = ' + n3
	d = 'D = ' + n4
	

	if n1_ >= 5.0:

		temp = '{:4}    {:5}      {:10}   {:5}  {:5}\n'.format( obsid,orbid,a,detid_a,pixid_a)
		file.write(temp)
	if n2_ >= 5.0:
		temp = '{:4}    {:5}      {:10}   {:5}  {:5}\n'.format( obsid,orbid,b,detid_b,pixid_b)
		file.write(temp)

	if n3_ >= 5.0:
		temp = '{:4}    {:5}      {:10}   {:5}  {:5}\n'.format( obsid,orbid,c,detid_c,pixid_c)
		file.write(temp)

	if n4_ >= 5.0:	
		temp = '{:4}    {:5}      {:10}   {:5}  {:5}\n'.format( obsid,orbid,d,detid_d,pixid_d)
		file.write(temp)

	# --------------to find data gap ---------------
	









print(' \n \n Done. Check all_noise_and_tele_error.dat for details.\n\n')

file.write('\n\nThanks,\nRahul\n')
file.flush()
file.close()





	
























