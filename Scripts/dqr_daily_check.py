#! /usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from urllib.request import urlopen,Request
from bs4 import BeautifulSoup
import re
import os
from astropy.time import Time
import datetime as dt
from humanfriendly import format_timespan


parser = argparse.ArgumentParser()
parser.add_argument("start", help="Enter starting row",type=int)
parser.add_argument('stop', help="Enter last row", type=int)
args = parser.parse_args()


print('\n####################################  Calculating noise ratio and checking telemetery error   ####################################\n')

file = open('noise_ratio.dat','w')
file.write("Hi,\n \n")
file.write('---------------------------------------------\n')
file.write('                  Noise Ratio   \n')
file.write('---------------------------------------------\n')
file.write('ObsID   OrbitID    Quadrant    DetID   PixID   \n')
file.write('---------------------------------------------\n') 


   
url = 'http://www.iucaa.in/~astrosat/czti_dqr/'
req = Request(url, headers={'User-Agent': 'Mozzilla/9.0'})
html = urlopen(req).read()
webpage = html.decode('utf-8')
s = BeautifulSoup(webpage, 'html.parser')


links = s.find_all('a',href = True, limit = args.stop, target="_blank")


order =[]
arranged_links = []
for i in range(len(links)):
	t = []
	if len(links[i].get('href')) == 54:
		temp = links[i].get('href')[26:30]
		t.append(int(temp))
		temp = links[i].get('href')[38:-11]
		t.append(int(temp))
		order.append(t)
		temp = links[i].get('href')
		arranged_links.append(temp)


sort_order = sorted(order)
sort_index = sorted(range(len(order)), key=lambda k: order[k],reverse = True)


final_links=[]
for i in sort_index:
	final_links.append(arranged_links[i])



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
	
def convert_time(time):
    try: 
        time_formatted = dt.datetime.strptime(time[:-3], '%Y-%m-%d  %H:%M:%S.%f')    
    except:
        time_formatted = dt.datetime.strptime(time[:-3], '%Y-%m-%d  %H-%M-%S.%f')
	
    time_formatted = dt.datetime.timestamp(time_formatted)
    
    return time_formatted

			
utc_all = []
all_time = []
all_id = []		
	
loop_start = len(final_links) 


for l in range(loop_start -1,0,-1):
	
	
	strip = final_links[l]
	#strip = links[l].get('href')
	link = 'https://www.iucaa.in/~astrosat/czti_dqr/' + strip
	html_ = urlopen(link)
	
	soup = BeautifulSoup(html_, 'lxml')
	
	info = soup.find_all('li')[4]
	obsid = info.text[-5:-1]
	

	orb = soup.find_all('title')[0]
	orbid = orb.text[-5:]
	

	
	ids = [obsid,orbid]
	all_id.append(ids)

	start_time = soup.find_all('li')[0].text[11:-1]
	start_time = start_time +' ' + soup.find_all('li')[1].text[11:-1]
	
	end_time = soup.find_all('li')[2].text[11:-1]
	end_time = end_time + ' ' + soup.find_all('li')[3].text[11:-1]
	
	utc_present = [start_time,end_time]
	
	start_time = convert_time(start_time)
	end_time = convert_time(end_time)
	
	present = [start_time,end_time]
	
	utc_all.append(utc_present)
	all_time.append(present)


	print('    Checking ObsID:', obsid ,	' OrbID:', orbid)

	
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
		print ('\n        Telemetry error is present in {} ,  orbit {}\n'.format(info.text[:-1],orbid))

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



print(' \n \n Done. Check noise_ratio.dat for details.\n\n')


print('####################################  Checking for data gap and ignore orbits   ################################################## \n')

flag = 0
all_miss = []
ignore_orbits = []

for i in range(1,len(all_time)):
	
	compare = all_time[i][0] - all_time[i-1][1]
		
	if (compare) > 0.0 :
		
		time = format_timespan(compare)
		if int(all_id[i][0]) == int(all_id[i-1][0]):
			string = "Data gap in observation {} between orbit {} and {} of ~ {}. \n".format(all_id[i][0], all_id[i-1][1],all_id[i][1],time)
			print(string) 
			flag = 1
			miss = [int(all_id[i-1][1]),int(all_id[i][1])]
			all_miss.append(miss)
		elif int(all_id[i][0]) != int(all_id[i-1][0]) and (compare) > 1800.0 :
			
			string = "Data gap between observation  {} and {} of ~ {}.\n".format(all_id[i][0], all_id[i-1][0],time)
			print(string)
			flag = 1
			miss = [int(all_id[i-1][1]),int(all_id[i][1])]
			all_miss.append(miss)
	
		#elif int(all_id[i][0]) != int(all_id[i-1][0]) and (compare) < 1800.0 :
		#	continue
	
	
	
	else :
		j = i + 1 
		while j < (len(all_time)):
			

			temp = all_time[j][0] - all_time[i-1][1]
			
			if temp < 0.0 :
				
				print("Orbit {} of observation {} can be ignored.\n".format(all_id[j-1][1],all_id[j][0]))
				ignore_orbits.append(all_id[j-1][1])
				j = j + 1
				flag = 1
			else :
				break




missing_orbit = []

for miss in all_miss:
	i = miss[0]
	
	while i < miss[1]:
		i = i + 1	
		missing_orbit.append(str(i))	
	missing_orbit.pop(-1)
	
string = ', '.join(missing_orbit)


if len(string) != 0:

	print ( "Missing orbits : ", string, '\n')

string = ', '.join(ignore_orbits)

if len(string) != 0:
	print ( "     Ignore orbits : ", string, '\n')

if flag == 0:
	print (' No data gap and ignore orbits. \n')


file.write('\n\nThanks,\nRahul\n')
file.flush()
file.close()





