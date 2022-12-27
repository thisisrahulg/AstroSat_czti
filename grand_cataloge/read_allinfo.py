#! /usr/bin/env python3
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--num", type = int, help = 'number of GRBs to select from all info')
args = parser.parse_args()

download_allinfo_cmd = 'sshpass -pC@dm1um rsync -zavh cztipoc@192.168.11.37:/data2/czti/special/GRB_ALL_info.dat .'
#download_allinfo_cmd = 'rsync -zavh cztipoc@192.168.11.37:/data2/czti/testarea/vipul/grb_all_info.txt GRB_ALL_info.dat'
os.system(download_allinfo_cmd)

grb_file = open("GRB_ALL_info.dat", "r")

lines = grb_file.readlines()
#lines =['GRB211105A 33009 299.53064 35.22754 81.93 -57.38 2021-11-05 04:35:20.20  373782921.20000005 -162.83717 90.50043 86.96948 28.34345 27.18598 632.00507 -5.19733 103.24031   VISIBLE  146.58634  109.51284 -167.57298  148.12539']

if args.num == None:
    num_rows = len(lines)
else:
    num_rows = args.num


all_info = []

for i in range(num_rows):
    all_info.append(lines[-i].split())

grb_file.close()

all_info = pd.DataFrame(all_info, columns = ["GRB_No","Orbi_No","RA_pnt","DEC_pnt","RA","DEC","Date","UTC_time","Astrosat_time","Rotation","Sun_angle","Moon_angle","Earth_angle","Elevation","alt","Earth_lat","Earth_long","Vis_ocu","Theta","Phi","Theta_x","Theta_y"])

all_info['space'] = ' '
all_info = all_info.sort_values(['Astrosat_time'], ascending = [True])

data_catalog = all_info[["GRB_No","space", "UTC_time", "Astrosat_time", "Orbi_No", "RA_pnt", "DEC_pnt", "Rotation", "Sun_angle", "Moon_angle", "Earth_angle", "alt", "Earth_lat", "Earth_long", "Theta","Phi", "Theta_x", "Theta_y", "RA", "DEC", "space", "Vis_ocu", "Date", "Elevation"]]

data_search = all_info[["GRB_No", "Date","UTC_time", "RA","DEC", "Vis_ocu","Theta","Phi", "Elevation", "Orbi_No"]]

data_catalog.to_excel("grb_all_sample.xlsx", index = False)
#data_search.to_csv("grb_search.txt", sep = ' ', index = False)
