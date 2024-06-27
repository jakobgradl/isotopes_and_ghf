
# AgeToIso.py

# Convert basal age to basal ice d18O composition
# Based on a linear correlation between Ngrip d18O record and Vostok CO2 record

# Need to provide the age of the basal ice (output of  NyePlusMelt.py)
# Make sure to update the path to the CO2 data

# ---------

# python packages
import numpy as np
import csv
import math as m
from copy import deepcopy
from netCDF4 import Dataset

# ISSM packages
from model import *

# main
def AgeToIso(LoadedModel, Age):
    # LoadedModel is your ISSM model; Age is the basal ice age from NyePlusMelt.py
    
    md = LoadedModel
    age = Age

    # load vostok co2 data
    co2Data = open('../Data/IsotopeData/Petit99_CO2_copy.txt','r')
    data = [i for i in csv.reader(co2Data,delimiter='\t')] # data has two columns: gas_age, ppm_co2
    data = data[1:] # get rid of header

    # create needed data arrays
    isocomp = np.zeros(md.mesh.numberofvertices)
    co2extrapol = np.zeros(md.mesh.numberofvertices)

    # match the basal ice age with the associated co2 value at each mesh vertex
    # the vostok record only has data at discrete gas-ages, the basal ice age is computed on a continuous scale
    # meaning: they don't line up
    # here, we find the two gas-ages in the vostok record that the ice age lies in between of
    # the co2 value associated with the ice is then computed from the co2 data of the two vostok data points by means of linear interpolation
    for i in np.arange(md.mesh.numberofvertices):
        for j in np.arange(len(data)):
            if age[i] < 2342 or age[i] > 414085: # that's the youngest/oldest gas-age in vostok record
                co2extrapol[i] = 280.0 # for ages younger than 2342 yrs, 280ppm is a good value; ages above 414085yrs aren't relevant in our study
                break
            elif age[i] < int(data[j][0]): # find vostok gas-age below vertex ice age
                continue
            elif age[i] >= int(data[j][0]): # find vostok gas-age above vertex ice age
                index = j
                dist = age[i] - int(data[index][0])
                gap = int(data[index+1][0]) - int(data[index][0])
                weight = dist / gap
                co2extrapol[i] = float(data[index][1]) + weight * (float(data[index+1][1]) - float(data[index][1])) # weighted mean/linear interpolation

    isocomp = co2extrapol * 0.08424944 - 58.729291
    # linear correlation between the vostok co2 record and the ngrip d18O record treated with a 500yr running mean
    isocomp = np.where(isocomp == -58.729291, np.nan, isocomp)
    
    return isocomp # d18O composition of the basal ice


