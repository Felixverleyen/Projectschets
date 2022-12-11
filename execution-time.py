import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

N, t = np.loadtxt('execution-times.txt', unpack = True)


#xdata = np.linspace(2,7,100)
#y=xdata**2
#plt.plot(xdata,y)

def func(x,a,b,c):
    return a*x**2 + b*x + c

par, var = curve_fit(func, N, t)
print(par)

plt.plot(N, func(N, *par), color='c', label = 'Expected execution time ~O(N²)')
plt.plot(N,t, color='b', label = 'Measured execution time')
plt.scatter(N,t, color='b')
plt.title('Mean duration time per integration time step in RK4 in function of the number of bodies')
plt.xlabel('# bodies')
plt.ylabel('Mean duration time [microseconds]')
plt.legend()
plt.show()




