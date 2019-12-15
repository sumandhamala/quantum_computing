"""Created on Fri Jul 13 18:19:18 2018 @author: s091d230
This code is to find the normalized rate of transition from one well of 1D CJJ
 qubit potential to other well using Montecarlo method. The barrier height of 
 the potential is fixed 0.29*2pi and the qubit flux is 0.5*2pi. The variable 
 quantity is exponent which fixed the Temperature(thermal noise term)
"""
# import the useful libraries for the simulation
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
#from lmfit import Model
#import seaborn as sea
import operator
import math 
pi = math.pi
##
dtau = (2*pi/ 150)   # time step
dtauinv = 1/dtau     # inverse of time
Phix = 0.5*2*pi      # External tilt knob
phixdc = 0.29*2*pi # expt value of cjj flux fixing the barrier
################ Langavin equation parameters ###########
beta0 = 2.43 
inv_bo = 1/beta0 
L = 100e-12  
K_B = 1.3807e-23 
U_0 = (2.068e-15)**2/((4*pi**2)* L) 
ej = ( U_0*beta0 )  
noiseAmp =( U_0*beta0 )/K_B  # express in unit of k (Ej and dU should be nearly equal)     
const2 = U_0/K_B 
const = beta0* np.cos( 0.5*phixdc) # potential const

''' Finding position of the particle in the left well '''  
## derivative of  1 dim. potential
def derv(x):
    return ( x - Phix ) + const* np.sin(x)  
## Make a starting guess 
x0 = np.array([1.6,3,4.65 ])  # for phix=0.5/ phixdc = 0.29
position = fsolve( derv, x0 )
## finding potential energy value
def potenergy(X,phix):
    return const2*( 0.5*( X- phix )**2  - const* np.cos(X) )
Uvalue = potenergy( position, Phix )
## finding inflectin point and threshold
def inflection(x):
    return 1 + const*np.cos(x)
X0 = [1, 4]
Inflection = fsolve( inflection, X0 ) 
##############################################################################
phi_L = position[0]
phi_R = position[2]
phi_saddle = position[1]
phi_plus = Inflection[1]   # right threshold
ddphi = phi_plus - phi_saddle 
phi_minus = phi_saddle - ddphi   # left threshold
DU0 =  Uvalue[1]- Uvalue[0]  #barrier height

G = 1.03           
  
noiseAmp = const2*beta0  #( U_0*beta0 )/k ; % express in unit of k (Ej and dU should be nearly equal)     
dtau2 = dtau**2 
C6 = inv_bo*dtau2
C44 = np.sqrt( (24*dtau**3)/noiseAmp)
C1 = (2-C6)-G*dtau
C2 =(G*dtau-1)
const1 = dtau2*np.cos( 0.5*phixdc ) ;
###############################################################################
    ## Changing parameters
Exp = np.arange(3,6,1)  #11:0.5:20;
jumpmax= 10
    ##
## initialize
downtime = np.zeros( (jumpmax, len(Exp) ) ) 
uptime = np.zeros( (jumpmax, len(Exp)) ) 
T = np.zeros( (len(Exp),) )
 
#############################################################################
for ii in range(len(Exp)):
    ##
    T[ii] = DU0/Exp[ii]
    C4 = C44*np.sqrt(np.multiply(G,T)) 
    phi1= phi_L; dphi= 0.00000001;
    phi2= phi1 + dphi;
    ##
for jump in range(0,jumpmax):
    count=0;
    randnumber = np.random.uniform(-0.5,0.5,1) 
    phi3 = C6*Phix + C1*phi2 + C2*phi1 - const1* np.sin(phi2) + C4[ii] *randnumber
    ## Initialized
    phi1 = phi2 #order matters
    phi2 = phi3
##
    while phi2 < phi_plus:
        ##
        count = count+1
        randnumber = np.random.uniform(-0.5,0.5,1)
        phi3 = C6*Phix + C1*phi2 + C2*phi1 - const1* np.sin(phi2) + C4[ii] *randnumber
        phi1 = phi2; ## order matters
        phi2 = phi3;
        ##
        downtime[jump,ii] = count*dtau    
        count=0
        ##
    while phi2 > phi_minus : 
        count = count+1
        randnumber = np.random.uniform(-0.5,0.5,1)
        phi3 = C6*Phix + C1*phi2 + C2*phi1 - const1* np.sin(phi2) + C4[ii] *randnumber
        phi1 = phi2  ## order matters
        phi2 = phi3 
        ##
        uptime[jump,ii] = count*dtau   
        ##















