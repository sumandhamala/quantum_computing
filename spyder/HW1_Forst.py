
# import essential libraries
import matplotlib.pyplot as plt
import numpy as np
pi= np.pi

time = np.linspace(0,100,1000)
y1 = 2* np.cos(2*pi*50*time)
y2 = 2* np.cos(2*pi*50*(time-0.0025))
y3 = 2* np.cos(2*pi*50*time- pi/4)



plt.figure(figsize=(8, 6), dpi=80)
plt.plot( time , y1, color="blue", linewidth=1.5, linestyle="-" )
plt.plot( time, y2, color="red", linewidth=1.5, linestyle="-" )
plt.plot( time, y3, color="green", linewidth=1.5, linestyle="-" )

#plt.legend()
plt.title('State',{'fontsize':20})
#plt.xlim([-2,2])
plt.xlabel('Time (s)',{'fontsize':20})
plt.ylabel('Amplitude',{'fontsize':20})
plt.tick_params(axis='both',labelsize=20)
plt.show()