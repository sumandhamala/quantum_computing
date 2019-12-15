""" Created on Mon Aug 13 14:15:20 2018@author: Suman Dhamala
 This code is to solve the Langevin equation for one dimensional
 flux qubit using RK4 method. This code is to solve any second order differential 
equation using RK4 method """
##################################################################################
# import the useful libraries for the simulation
import numpy as np
import matplotlib.pyplot as plt
import time
import math 
pi = math.pi

#################################################################
dt = 2*pi/150.0     # Time step
G  = 1.03           # dimensionless damping parameter (1/W_pRC)
beta = 2.43         #
inv_beta = 1/beta;
phixdc = 0.29*2*pi  #
coeff = 1/6.0 
coeff1 = np.cos(0.5*phixdc) 
#################################################################
nu = 10**-5.8          # ramp rate
m = round ((0.5*2*pi - 0.465*2*pi)/(nu*dt) )  # no. of data points for slow ramp
mm = round(m/1.5)       #flat part
mmm = 70                #fast ramp 
feedback = ( m + mm ) 
##Creating External flux 
slow_ramp =  np.linspace(0.465*2*pi, 0.5*2*pi, m  )
flat_part = np.full( (mm, ), 0.5*2*pi, dtype=float)
fast_ramp = np.linspace(0.5*2*pi , 0.535*2*pi,mmm)
endpoint = np.array([0.535*2*pi+0.00001] )
fast_ramp1 = np.linspace(0.5*2*pi , 0.465*2*pi,mmm)
endpoint1 = np.array([0.465*2*pi-0.00001] )
##
phix = np.concatenate((slow_ramp, flat_part, fast_ramp, endpoint))
#phix1 = np.concatenate((slow_ramp, flat_part, fast_ramp1, endpoint1))

## uniform distributed random number to create noise term  
randnumber = np.random.uniform(-0.5,0.5,len(phix)) 
###############################################################
phi = np.zeros((len(phix),))
v = np.zeros((len(phix),))
t = np.zeros((len(phix),))
phi[0]= 1.4442            # initial position
v[0] = 0              #initial velocity

###############################################################
def func( phi, v, phixx):
    return -G*v + inv_beta*(phixx-phi) - coeff1*np.sin(phi)

################## 4th order Rungakutta method #################
t1 = time.time() #times the computation
count = 0   
while phix[count]  <  phix[-1]: #phix[round(feedback+1 )]:
    a1 = dt * func(phi[count], v[count] , phix[count])
    b1 = dt * v[count]
    a2 = dt * func( phi[count]+0.5*b1 , v[count]+0.5*a1 , phix[count]+0.5*a1)
    b2 = dt * (v[count]+0.5*a1)
    a3 = dt * func( phi[count]+0.5*b2 , v[count]+0.5*a2 , phix[count]+0.5*a2)
    b3 = dt * (v[count]+0.5*a2)
    a4 = dt * func( phi[count]+b3 , v[count]+a3 , phix[count]+a3)
    b4 = dt * (v[count]+a3)
    
    phi[count+1] = phi[count]+coeff*(b1+2*(b2+b3)+b4)+0.003*randnumber[count]
    v[count+1] = v[count]+coeff*(a1+2*(a2+a3)+a4)+0.003*randnumber[count]
    #t[count+1]= t[count]+dt
    count += 1
t2 = time.time()
print ('computation takes ',t2-t1,' seconds.')

plt.plot( phi[1:-1] )
#plt.plot( v[1:-1] )









#t1 = time.time() #times the computation
#count = 0   
#while phix[count]  <  phix[round(feedback+1 )]:
#    k1x = v[count]
#    k2x = v[count] +dt * k1x/2
#    k3x = v[count] +dt * k2x/2
#    k4x = v[count] +dt * k3x
#    
#    k1v = func(phi[count], v[count] , phix[count])
#    k2v = func(phi[count]+dt*k1x/2, v[count]+dt*k1v/2 , phix[count])
#    k3v = func(phi[count]+dt*k2x/2, v[count]+dt*k2v/2 , phix[count])
#    k4v = func(phi[count]+dt*k2x, v[count]+dt*k2v , phix[count])
#    
#    phi[count+1] = phi[count]+ dt*coeff*( k1x+2*(k2x+k3x)+k4x )+0.003*randnumber[count] 
#    v[count+1]   = v[count]+ dt*coeff*( k1v+2*(k2v+k3v)+k4v )
#    count += 1
#    timee[count+1]= timee[count]+dt
#t2 = time.time()
#print ('computation takes ',t2-t1,' seconds.')
#
#
#
#
#t1 = time.time() #times the computation
#count = 0   
#while phix[count]  <  phix[round(feedback+1 )]:
#    k1x = dt * v[count]
#    k2x = dt * ( v[count] + k1x/2  )
#    k3x = dt * ( v[count] + k2x/2 )
#    k4x = dt * ( v[count] + k3x )   
#    
#    k1v = dt *func(phi[count], v[count] , phix[count])
#    k2v = dt *func(phi[count]+ k1x/2, v[count]+ k1v/2 , phix[count])
#    k3v = dt *func(phi[count]+ k2x/2, v[count]+ k2v/2 , phix[count])
#    k4v = dt *func(phi[count]+ k2x, v[count]+ k2v , phix[count])
#    
#    phi[count+1] = phi[count]+ coeff*( k1x+2*(k2x+k3x)+k4x )+0.003*randnumber[count] 
#    v[count+1]   = v[count]+ coeff*( k1v+2*(k2v+k3v)+k4v )
#    count += 1
#t2 = time.time()
#print ('computation takes ',t2-t1,' seconds.')
