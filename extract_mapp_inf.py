#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extracts all Log.final.out infos from STAR count mode output

@author: johaGL 2022
"""
import sys
import os
import pandas as pd
import yaml

configfile = sys.argv[1]
outdir = sys.argv[2]

def open_config_file_snake_version(confifile):
    try:
        with open(confifile, "r") as f:
            confidic = yaml.load(f, Loader=yaml.Loader)            
    except yaml.YAMLError as yam_err:
        print(yam_err)
        print("\nimpossible to read configuration file")
        confidic = None
    except Exception as e:
        print(e)
        print("\nimpossible to read configuration file")
        confidic = None

    return confidic
    
confidic = open_config_file_snake_version(configfile)
# --

print("mapping final reports (mapping_infos.csv), into directory :", outdir)

locationsams = confidic["mappeddir"]
samplesfile = confidic["metadata"]

df = pd.read_csv(samplesfile)

mysamples = list(df['sample'])


myfi_ =  [ locationsams+samp+"/Log.final.out" for samp in mysamples] 


dem =  {}
        
for k in myfi_:
    with open(k,'r') as f:
        fli = f.readlines()
    for li in fli:
        elef = li.replace('\n','').split("|")
        if len(elef) > 1:
            mykey = elef[0].replace("   ","").strip()
            myval = elef[1].replace("\t","").strip()
            if myval != "":
                if mykey not in dem:
                    dem[mykey]= []
                    dem[mykey].append(myval.replace("%",""))
                else:
                    dem[mykey].append(myval.replace("%",""))
                

dem2 = {}
for key in dem.keys():
    if key  == "Number of input reads" or key == "Average input read length":
        dem2[key] = dem[key]
    elif "%" in key:
        dem2[key] = dem[key]
            

idf = pd.DataFrame.from_dict(data=dem, orient='index') 
idf2 = pd.DataFrame.from_dict(data=dem2, orient="index")
idf.columns = mysamples #replace by mysamples
idf2.columns = mysamples
idf.to_csv(f"{outdir}mapping_infos_complete.csv")
idf2.to_csv(f"{outdir}mapping_infos_reduced.csv")

