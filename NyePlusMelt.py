
# NyePlusMelt.py

# Calculate the age of the basal ice with the Nye+melt age model
# You will need to provide an estimate of basal meltrate

#-----------


# python packages
import numpy as np
from importlib import reload
from copy import deepcopy

# ISSM packeges
from model import *

# main
def nyeplusmelt(LoadedModel, Melt): # LoadedModel is your ISSM model
    md = deepcopy(LoadedModel)
    M = deepcopy(Melt)

    H = md.geometry.thickness
    z = md.mesh.z - md.geometry.base
    A = md.smb.mass_balance

    enum = (A-M)*H + H*M
    denom = (A-M)*z + H*M
    denom = np.where(denom==0, np.nan, denom)
    
    a = (H / (A-M)) * np.log(enum/denom)

    return a
