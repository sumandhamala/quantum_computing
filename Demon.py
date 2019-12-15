''' Coded by suman Dhamala
#life time experiment for demon '''
####################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
from datetime import *
#####################################################
for i in range(1,21):
    fname='demon_%d_Vin(0.6500V)_Vsq(1.2017V)_Vrf(-3.0000V)_Vcjj(-5.1000V).dat'%(i)
data = pd.read_csv(fname , sep='\t' , header = None , skiprows = 0 , comment = '#')
Vqubit = data.iloc[:, [0]].values 
Vout = data.iloc[:, 1:].values
Vout_reshape = Vout.reshape(Vout.shape[0]*10,10)



#Vout = Vout.ravel(2, C)
#Vout = np.concatenate(d2,10)
#A = np.squeeze(np.asarray(d2))
#d2.shape
#Vout = d2.ravel()