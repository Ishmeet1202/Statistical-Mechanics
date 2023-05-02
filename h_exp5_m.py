import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.stats import linregress

k = 1.38e-23
h = 6.624e-34
c = 3e8


def planck(x):
    o = x**3/(np.exp(x)-1)
    return o

planck = np.vectorize(planck)

ino = 1e-2 ; ut = 15
x = np.linspace(ino,ut,100)

def mean(tol,pt):
    ar = integrate.simps(planck(x),x)
    x_n = np.linspace(ut/8,ut/2,pt)
    for i in range(len(x_n)):
        x_nn = np.linspace(ino,x_n[i],100)
        ar_n = integrate.simps(planck(x_nn),x_nn)
        if abs(ar_n - ar/2 ) < tol :
           return x_n[i]
           break 
             
lm_m = mean(0.01,500) 
print('x_median',lm_m)

b = h*c/(k*lm_m)
print('wein constant',b)
    
# stefan

e = lambda t : k*T
l_p = lambda t : h*c/e(t)

I_p = integrate.quad(planck,1e-15,20)

print(I_p[0],np.pi**4/15)

def u_d(T):
    return (np.pi**4/15)*(8*np.pi*e(T)**4/(h*c)**3)

T = np.arange(100,10500,500)

F = lambda T : (c/4)*u_d(T)
    
res = linregress(np.log(T),np.log(F(T)))

print('slope',res[0],'\n','intercept',res[1])

print('stefan constant', np.exp(res[1]) )

fig,ax = plt.subplots()
ax.plot(x,planck(x))
ax.set_title('Planck law of radiation')
ax.set_xlabel('x')
ax.set_ylabel('energy spectral density')
plt.show()

fig,ax = plt.subplots()
ax.plot(T,F(T))
ax.set_title('Radiant flux vs T')
ax.set_xlabel('T')
ax.set_ylabel('Radiant flux')
plt.show()

fig,ax = plt.subplots()
ax.plot(np.log(T),np.log(F(T)),'--o')
ax.plot()
ax.set_title('log(F)')
ax.set_xlabel('log(T)')
ax.set_xlim([-1,9])
ax.set_ylim([-20,20])
plt.grid()
ax.set_ylabel('regression')
plt.show()



