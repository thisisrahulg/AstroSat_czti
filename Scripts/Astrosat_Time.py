#!/usr/bin/env python2.7

"""
    Astrosat_Time.py
    Aim: This is a time conversion routine for Astrosat CZTI Instrument, It can be used to convert Astrosat Seconds (since epoch, 2010.0 UTC) to ISO 8601 date,
    Calender Date, Year and DayNumber, Julian Day and Modified Julian Day formats and Vice Versa.
    
    Type Astrosat_Time.py -h to get help

"""

#------------------------------------------------------------------------

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.time import Time
from datetime import datetime


##### Defining default values and constants

Epoch_Time=157766400; # Total numbers of seconds from 2010.00.00.00.00 till 2015.00.00.00.00
Epoch_Date=Time('2010-01-01 00:00:00', scale='utc')

def convertISODatetosec(ISO_Date):
        Astrosat_second=ISO_Date-Epoch_Date
        Astrosat_second=Astrosat_second.sec-3;
        Astrosat_iso=ISO_Date.iso
        Astrosat_yday=ISO_Date.yday
        Astrosat_jd=ISO_Date.jd
        Astrosat_mjd=ISO_Date.mjd
        return Astrosat_second + 1.0

#print convertISODatetoSec(Current_Date)        

def convertYDAYtosec(YDAY):
        Epoch_Date_yday=Time(Epoch_Date,format='yday')
        Astrosat_second=(YDAY-Epoch_Date_yday).sec-3;
        Astrosat_yday=YDAY.yday
        Astrosat_iso=YDAY.iso
        Astrosat_jd=YDAY.jd
        Astrosat_mjd=YDAY.mjd
        return Astrosat_yday,Astrosat_second,Astrosat_iso,Astrosat_jd,Astrosat_mjd
#print convertYDAYtosec(Current_yday)


def convertJDtosec(JD):     
        Epoch_Date_jd=Time(Epoch_Date,format='jd')
        Astrosat_second=(JD-Epoch_Date_jd).sec-3;
        Astrosat_jd=JD.jd
        Astrosat_iso=JD.iso
        Astrosat_yday=JD.yday
        Astrosat_mjd=JD.mjd
        return Astrosat_jd, Astrosat_second, Astrosat_iso, Astrosat_yday, Astrosat_mjd
#print convertJDtosec(Current_jd)


def convertMJDtosec(MJD):
       Epoch_Date_mjd=Time(Epoch_Date,format='mjd')
       Astrosat_second=(MJD-Epoch_Date_mjd).sec-3;
       Astrosat_mjd=MJD.mjd
       Astrosat_iso=MJD.iso
       Astrosat_yday=MJD.yday
       Astrosat_jd=MJD.jd
       return Astrosat_mjd, Astrosat_second,Astrosat_iso,Astrosat_yday,Astrosat_jd
#print convertMJDtosec(Current_mjd)


def convertAStoAll(Astrosat_second):
	ND=int((Astrosat_second-Epoch_Time)//86400); # No of days  
	NS=(Astrosat_second-Epoch_Time)%86400;
	NH=int(NS//3600);  # No of hours 
	NHS=NS%3600;
	NM=int(NHS//60); # No of minutes 
	#print NM
	NOS=float(NHS%60); # No of seconds  with fraction 
	#print NOS
	NOSI=str(NOS)[0]
	NOSF=str(NOS)[1:]
	year=2015;
	dayr=ND;
	while (dayr>0):
		if((year%4)==0):
			NDY=366
		else:
			NDY=365
		dayr=dayr-NDY;
		if (dayr>0):
			year=year+1;
		else: 
			DaY=dayr+NDY+1;
	
	if NH==0 and NM==0 and NOS==0.0:
		Astrosat_Yr=str(year)+str(':')+str(DaY)+str(':')+str(NH)+str(NH)+str(':')+str(NM)+str(NM)+str(':')+NOSI+NOSI+NOSF;
	elif NH==0:
		Astrosat_Yr=str(year)+str(':')+str(DaY)+str(':')+str(NH)+str(NH)+str(':')+str(NM)+str(':')+NOSI+NOSF;
	elif NM==0:
		Astrosat_Yr=str(year)+str(':')+str(DaY)+str(':')+str(NH)+str(':')+str(NM)+str(NM)+str(':')+NOSI+NOSF;
	elif NOS==0.0:		
		Astrosat_Yr=str(year)+str(':')+str(DaY)+str(':')+str(NH)+str(':')+str(NM)+str(':')+NOSI+NOSI+NOSF;
	elif NH==0 and NM==0:
		Astrosat_Yr=str(year)+str(':')+str(DaY)+str(':')+str(NH)+str(NH)+str(':')+str(NM)+str(NM)+str(':')+NOSI+NOSF;
	elif NM==0 and NOS==0.0:
		Astrosat_Yr=str(year)+str(':')+str(DaY)+str(':')+str(NH)+str(':')+str(NM)+str(NM)+str(':')+NOSI+NOSI+NOSF;
	elif NH==0 and NOS==0.0:
		Astrosat_Yr=str(year)+str(':')+str(DaY)+str(':')+str(NH)+str(NH)+str(':')+str(NM)+str(':')+NOSI+NOSI+NOSF;	   
	else: 
		Astrosat_Yr=str(year)+str(':')+str(DaY)+str(':')+str(NH)+str(':')+str(NM)+str(':')+NOSI+NOSF;

	#print (Astrosat_Yr)
	Astrosat_yday_obj=Time(Astrosat_Yr, scale='utc')
	Astrosat_yday=Astrosat_yday_obj.yday
	Astrosat_iso=Astrosat_yday_obj.iso
	Astrosat_jd=Astrosat_yday_obj.jd
	Astrosat_mjd=Astrosat_yday_obj.mjd
	return Astrosat_second,Astrosat_yday,Astrosat_iso,Astrosat_jd,Astrosat_mjd       
#print convertAStoAll(206021562.0)

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
group.add_argument("-s",     "--Astrosat_second", help= "Enter ASTROSAT seconds since 2010.0 UTC (decimal) to covert to other formats ", type=float)
group.add_argument("-iso",   "--Astrosat_iso",    help= "Enter Astrosat Time in ISO 8601 date (yyyy-MM-dd hh:mm:ss) to covert to other formats ", type=str)
group.add_argument("-yday",  "--Astrosat_yday",   help= "Enter Astrosat Time in Year and day number (yyyy:ddd:hh:mm:ss) to covert to other formats ", type=str)
group.add_argument("-jd",    "--Astrosat_jd",     help= "Enter Astrosat Time in Julian Day (ddddddd.ddd...) to covert to other formats ",  type=float)
group.add_argument("-mjd",   "--Astrosat_mjd",    help= "Astrosat Time in Modified Julian Day (ddddd.ddd...) to covert to other formats  ", type=float)
args = parser.parse_args()

if args.Astrosat_second:
        print (convertAStoAll(args.Astrosat_second))

if args.Astrosat_iso:
       Current_iso=Time(args.Astrosat_iso, scale='utc')
       print (convertISODatetosec(Current_iso))

if args.Astrosat_yday:
       Current_yday=Time(args.Astrosat_yday, scale='utc')
       print (convertYDAYtosec(Current_yday))

if args.Astrosat_jd:
       Current_jd=Time(args.Astrosat_jd, format='jd', scale='utc')
       print (convertJDtosec(Current_jd))

if args.Astrosat_mjd:
       Current_mjd=Time(args.Astrosat_mjd,format='mjd', scale='utc')
       print (convertMJDtosec(Current_mjd))	












