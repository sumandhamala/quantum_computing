''' This code is to solve the Langevin equation for one dimensional flux qubit.
The maxwell demon experiment is realised through the monte Carlo simulation method.
In order to digitized the position of the particle in a 1D potnetial the 
Double Gaussian method is used '''
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
##
nu = 10**-5.5  # ramp rate
m = round ((0.5*2*pi - 0.465*2*pi)/(nu*dtau) )  # no. of data points for slow ramp
mm = round(m/1.5)   #flat part
mmm = 70  #fast ramp 
feedback = ( m + mm ) 
##Creating External flux 
slow_ramp =  np.linspace(0.465*2*pi, 0.5*2*pi, m  )
flat_part = np.full( (mm, ), 0.5*2*pi, dtype=float)
fast_ramp = np.linspace(0.5*2*pi , 0.535*2*pi,mmm)
endpoint = np.array([0.535*2*pi+0.00001] )
fast_ramp1 = np.linspace(0.5*2*pi , 0.465*2*pi,mmm)
endpoint1 = np.array([0.465*2*pi-0.00001] )
##
Phi_x = np.concatenate((slow_ramp, flat_part, fast_ramp, endpoint))
Phi_x1 = np.concatenate((slow_ramp, flat_part, fast_ramp1, endpoint1))
## creating Epsilon
Ep_slowramp = np.linspace(-12.0147, 0, m)
Ep_flatpart = np.zeros(mm,)
Ep_fastramp = np.linspace(0, 12.0147, mmm )
Ep_fastramp1 = np.linspace(0, -12.0147, mmm )
Ep_endpoint = np.array([12.0147+0.00001] )
Ep_endpoint1 = np.array([-12.0147-0.00001] )
##
Epsilon = np.concatenate( (Ep_slowramp, Ep_flatpart, Ep_fastramp, Ep_endpoint  ) )
Epsilon1 = np.concatenate( (Ep_slowramp, Ep_flatpart, Ep_fastramp1, Ep_endpoint1 ) )
##
phixdc =   0.29*2*pi  # expt. value of cjj flux
Totaltau = (m + mm + mmm +1)*dtau  #normalised total time
aka = nu*dtau 
## convert to real time
TotalT = np.arange(0,Totaltau,dtau)
Trealtime = ( TotalT )/10**2  # converted to nanosecond
slowramptime = (m)*dtau/10**2  
##
################ Langavin equation parameters ###########
beta0 = 2.43 
inv_bo = 1/beta0 
L = 100e-12  
K_B = 1.3807e-23 
U_0 = (2.068e-15)**2/((4*pi**2)* L) 
ej = ( U_0*beta0 )  
noiseAmp =( U_0*beta0 )/K_B  # express in unit of k (Ej and dU should be nearly equal)     
const2 = U_0/(K_B*ej) 
##
dtau2 = dtau**2 
C6 = inv_bo*dtau2
G = 1.03        
gamma = 18.64    
#Exp = 27.12  # for phicjj= 0, T = 4.2K 
Exp = 4.84  # for phicjj= 0.29, T = 4.2K 
##
D6 = C6*4*gamma;
C44 = np.sqrt( (24*dtau**3)/noiseAmp)
C1 = (2-C6)-G*dtau
C2 =(G*dtau-1)
## 
#DU0 = Uvalue(2)- Uvalue(1);  %% barrier height when phix 0.5
DU0 = 20.332 #for phicjj= 0.29, phix=0.5
#DU0 = 113.91 #for phicjj= 0, phix=0.5
T = DU0/Exp
C4 = C44*np.sqrt(G*T) 
#EpsilonT = Epsilon/T;
const4 = G*dtauinv*noiseAmp 
const1 = dtau2*np.cos( 0.5*phixdc ) 
const  = beta0* np.cos( 0.5*phixdc )
const3 = 1/T  
## uniform distributed random number to create noise term  
randnumber = np.random.uniform(-0.5,0.5,len(Phi_x)) 
## Initialization of array
phi = np.zeros((len(Phi_x),))
Work = np.zeros((len(Phi_x),))
dizWork = np.zeros((len(Phi_x),))
Qheat = np.zeros((len(Phi_x),))
##
rows = 2
NetWORK = np.zeros((rows, ))
state = np.zeros((len(Phi_x), ))
#RESULT = np.zeros((rows,3)).reshape(2,3) 
#RESULT = np.empty((rows, 3))
RESULT = []
NETWORK = []
## 
''' Finding position of the particle in the left well '''  
## derivative of  1 dim. potential
def derv(x):
    return ( x - Phi_x[0] ) + const* np.sin(x)  
## Make a starting guess 
x0 = np.array([1.4,3.66,4.26] ) #for phix = 0.47
## x0 = np.array([1.6,3,4.65 ])  # for phix=0.5/ phixdc = 0.2
position = fsolve( derv, x0 )
## finding potential energy value
def potenergy(X,phix):
    return const2*( 0.5*( X- phix )**2  - const* np.cos(X) )
## finding inflectin point and threshold
#def inflection(x):
#    return 1 + const*np.cos(x)
#X0 = [1, 4]
#Inflection = fsolve( inflection, X0 ) 
########################################################################
for i in  range(0,rows):
## initail values of the position of the particle in the potential
    phi[0] = position[0]
    phi[1] = phi[0] + aka
    Work[0] = potenergy( phi[0], Phi_x[0]) - potenergy( phi[1], Phi_x[1])
    count = 1
##
    while Phi_x[count]  <  Phi_x[round(feedback+1 )]:
     ## loop to find position of the particle   
     phi[count+1] = C6*Phi_x[count] + C1*phi[count]+ C2*phi[count-1]- const1* np.sin(phi[count])+ C4*randnumber[count] 
     ## work done by the engine
     Work[count] = potenergy( phi[count], Phi_x[count]) - potenergy( phi[count+1], Phi_x[count+1])
     count = count+1
## feedback process has two fast ramp part
    if phi[count] >  position[1]:
            while Phi_x[count]  < endpoint:
                ## loop to find position of the particle
                phi[count+1] = C6*Phi_x[count] + C1*phi[count-1]+ C2*phi[count]- const1* np.sin(phi[count])+C4*randnumber[count] 
                ## work done by the engine
                Work[count] = potenergy( phi[count], Phi_x[count]) - potenergy( phi[count+1], Phi_x[count+1])
                count = count+1
    else:
            while Phi_x1[count]  > endpoint1:
             ## loop to find position of the particle   
                phi[count+1] = C6*Phi_x1[count] + C1*phi[count]+ C2*phi[count-1]- const1* np.sin(phi[count])+C4*randnumber[count] 
                ## work done by the engine
                Work[count] = potenergy( phi[count], Phi_x1[count]) - potenergy( phi[count+1], Phi_x1[count+1])
                count = count+1
## finding the net workdone by the engine
cumsumWORK = Work.cumsum()
Network = -const3*ej*cumsumWORK[-1]
NETWORK.append(Network)
NW = np.vstack(NETWORK)
############### TO CALCULATE work and heat after digitization###############################
###### Gaussian method for digitization #########
## Smoothing and mapping the position to [0,1]
#phi = smooth (smooth(phi(1:end)) ) ;
mapphi = ( phi - phi.min() )/( phi.max()- phi.min() )
## plotting histogram without showing figure 
counts, bin_edges = np.histogram( mapphi, bins=120, density=True)
bincenter = bin_edges[1:len(bin_edges)] + np.diff(bin_edges) / 2
xdata = bincenter 
ydata = counts 
## figure out initial guess by finding amplitude and mean from the data
index0, Amp0 = max(enumerate(ydata[:int(len(ydata)/2)]), key=operator.itemgetter(1))
mu0 = xdata[index0]
index1, Amp1 = max(enumerate(ydata[int(len(ydata)/2):]), key=operator.itemgetter(1))
mu1 = xdata[index1+int(len(ydata)/2)]
initialguess = [ Amp0, mu0, 0.05, Amp1 , mu1, 0.05] 
##
def gauss(x,A,mu,sigma):
    return A*np.exp(-0.5*( (x-mu)/sigma )**2 )
def bimodal(x,Amp0,mu0,sigma0,Amp1,mu1,sigma1):
    return gauss(x,Amp0,mu0,sigma0) + gauss(x,Amp1,mu1,sigma1)
## aNOTHER METHOD OF FITTING LMFIT
#Dgaussianmodel = Model(bimodal)
#pars = Dgaussianmodel.guess(ydata, xdata=xdata)
#paras = Dgaussianmodel.make_params(6, 0.19, 0.05, 1, 0.8, 0.05)
#result = Dgaussianmodel.fit( ydata, xdata=xdata )
##
coeff,_ = curve_fit( bimodal, xdata, ydata, initialguess )
#sigma = np.sqrt(np.diag(cov))
## plotting for checking the gaussian distribution
 ## check from histogram
#plt.figure(i)
#plt.hist(mapphi, density=True, bins=120, alpha=0.9)
#plt.plot(xdata,bimodal(xdata,*coeff),color='red',lw=2,label='model')
#plt.legend('Double gaussian fit')
#print(coeff,'\n',sigma) 
 ## threshold voltage and finding root of quadratic equation
aa =  ( coeff[5]**2- coeff[2]**2  );
bb =  2*( -coeff[1] * coeff[5]**2 + coeff[4] * coeff[2]**2 ) ;
cc = - coeff[4]**2 * coeff[2]**2 + coeff[1]**2 * coeff[5]**2 - 2* coeff[2]**2 * coeff[5]**2  * np.log(coeff[0]/coeff[3])
p = [aa, bb, cc]
r = np.roots(p) 
vthh = r[0] 
vth = r[1] 
##
## convert position to digitized value state[-1,1]
for jj in range(0,len(mapphi)):
    if mapphi[jj] > vth:
        state [jj] = 1
    else:
        state[jj] = -1 
##
## Choose epsilon accord. to feedbackpoint
if state[-1] == 1:
    E = Epsilon 
else:
    E = Epsilon1 
############################calculate heat and work############################
for ii in range(0,len(state)):
    if state[ii]== state[ii-1]:
            dizWork[ii]= dizWork[ii-1]+ 0.5*state[ii]*(E[ii]-E[ii-1])
            Qheat[ii]= Qheat[ii-1]
    else:
            dizWork[ii]= dizWork[ii-1]
            Qheat[ii]= Qheat[ii-1] - 0.5*( state[ii]-state[ii-1] )*( E[ii]+E[ii-1] )*0.5   
#################################################### #############################          
## saving the results ##########
RESULTt = [ dizWork[-1], Qheat[-1] ] 
RESULT.append(RESULTt)
finalresult = np.vstack(RESULT)
## initialized
count=0
#############################################################################
plt.figure(figsize=(8, 6), dpi=80)
plt.plot( Trealtime,mapphi, color="blue", linewidth=1.0, linestyle="-" )
plt.plot( Trealtime, state, color="red", linewidth=1.0, linestyle="-" )
#plt.legend()
plt.title('State',{'fontsize':20})
plt.xlim([0,1200])
plt.xlabel('Time ($\mu$s)',{'fontsize':20})
plt.ylabel('state',{'fontsize':20})
plt.tick_params(axis='both',labelsize=20)
plt.show()