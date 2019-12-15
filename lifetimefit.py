import numpy as np
from scipy.optimize import curve_fit

def f(t,N,t0,y0,tau):

	return y0+N*np.exp(-(t-t0)/tau)


w = np.loadtxt("suman.txt")
t = w[:, 0]
d1 = w[:,1]
d = d1+1
sig = np.sqrt(d)

start = (1600000, 0,0, 0.005)




popt, pcov = curve_fit(f,t,d, sigma = sig,p0 = start,absolute_sigma=True)

print('popt=',popt)
print('pocv=',pcov)
  
Nexp = f(t, *popt)
r = d - Nexp
chisq = np.sum((r/sig)**2)
   
print("chisq =",chisq)

import matplotlib.pyplot as plt
plt.errorbar(t, d, yerr=sig, fmt = 'o', label='"data"')
plt.plot(t, Nexp, label='fit')
plt.legend()
plt.savefig('suman.pdf')

plt.show()

 
