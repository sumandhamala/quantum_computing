"""
Random signal and noise, HW 8/23/2018, Author: Suman Dhamala
"""
# import libraries
import math
import numpy as np
import matplotlib.pyplot as plt
pi= np.pi
###############################################################################
# plotting the rectangular function having width 0.5 and center at 3
def rect(t):
    return np.where(abs(t)<=0.5, 1, 0)

time = np.linspace(-5, 5, 1000)
x= np.sqrt(2)* ( rect( (time+3)/0.5 )- rect( (time-3)/0.5) )
#xy= rect(time)
absx = abs(x*x)

plt.figure(figsize=(8, 6), dpi=80)
plt.plot( time , x, color="blue", linewidth=1.5, linestyle="-" )
#plt.legend()
plt.title('Combination of rectangular signal',{'fontsize':20})
#plt.xlim([-2,2])
plt.xlabel('t(s)',{'fontsize':20})
plt.ylabel('x(t)',{'fontsize':20})
plt.tick_params(axis='both',labelsize=20)
plt.show()
###############################################################################
# plotting fourier transfer of rect signal in Decible unit
f = np.linspace(-5, 5, 10000)
xff = 2*(np.sqrt(2))*(0.5)* (np.sinc(pi*f*0.5) )* (np.sin(2*pi*f*3) )
xf = 20*np.log(abs(xff))
plt.plot(f , xf, color="red", linewidth=2 , linestyle="-")
#plt.legend()
plt.title('FT of rect signal in Decible dB',{'fontsize':15})
#plt.xlim([-2,2])
plt.xlabel('f(Hz)',{'fontsize':15})
plt.ylabel('20Log(abs(X(f)))',{'fontsize':15})
plt.tick_params(axis='both',labelsize=15)
plt.show()
###############################################################################
# plotting complex numbers in a plot 
L=[ 0.94 + 0.34j, -0.34 + 0.94j, -0.94 -0.34j, 0.34+0.94j]
X = [x.real for x in L]
Y = [x.imag for x in L]
#take each number to draw from center
for xx in range(len(X)):
        plt.plot([0, X[xx]],[0, Y[xx]],'ro-',label='python')
# draw x and y line         
plt.plot([-1.1, 1.1], [0, 0], 'k-', lw=1)
plt.plot([0, 0], [-1, 1], 'k-', lw=1)
plt.xlim((-1.1,1.1))
plt.ylim((-1,1))
#plt.title('',{'fontsize':15})
plt.ylabel('Imaginary',{'fontsize':15})
plt.xlabel('Real',{'fontsize':15})
plt.tick_params(axis='both',labelsize=15)
plt.grid()
plt.show()
###############################################################################
h = 10*np.exp(-time/4)
plt.plot( time , h, color="blue", linewidth=1.5, linestyle="-" )
#plt.legend()
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
plt.title('time respone of H2(f)'.translate(SUB),{'fontsize':20})
#plt.xlim([-2,2])
plt.xlabel('t(s)',{'fontsize':20})
plt.ylabel('h2(t)'.translate(SUB),{'fontsize':20})
plt.tick_params(axis='both',labelsize=20)
plt.show()
###############################################################################
mag = ( (10/4)*np.sqrt(1 + (2*pi*f)**2 ) )/( (1/4)**2 + (2*pi*f)**2 )
magnitude = 20*np.log(abs(mag))

angle = math.atan(-2*pi*f)


plt.plot( f , magnitude, color="blue", linewidth=1.5, linestyle="-" )
#plt.title('time respone of H2(f)'.translate(SUB),{'fontsize':20})
#plt.xlim([-2,2])
plt.xlabel('f(Hz)',{'fontsize':20})
plt.ylabel('Phase',{'fontsize':15})
plt.tick_params(axis='both',labelsize=20)
plt.show()




