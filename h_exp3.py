import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as intn

td_c = 343
td_e = 230
r = td_e/td_c

def sp_einstein(x):
    
    return ((1/x)**2)*(np.exp(1/x)/(np.exp(1/x)-1)**2)

sp_einstein = np.vectorize(sp_einstein)

def debye(x):
    
    u1 = -3*x/(np.exp(x)-1)
    u2 = 12/x**3
    ine = intn.quad(lambda x : (x**3/(np.exp(x)-1)),0,x)
    return u1 + u2*ine[0]  

debye = np.vectorize(debye)

var = np.linspace(1e-8,2,50)

def dos_d(v):
    if v <= 1/r:
       return v**2
    else:
       return 0

def dos_e(v):
    if abs(v-1) <= 0.02:
       return 1
    else:
       return 0

dos_d = np.vectorize(dos_d)
dos_e = np.vectorize(dos_e)

#plotting DOS
v = np.linspace(0,2,100)
print(v)


fig,ax = plt.subplots()
ax.plot(v,dos_e(v),label = 'Einstein theory')
ax.plot(v,dos_d(v)/dos_d(1/r),'-o',label = 'Debye theory')
ax.set_ylim([0,1.1])
ax.set_xlabel('f/fe')
plt.legend(loc = 4)


#plotting
fig,ax = plt.subplots()
ax.plot(var,sp_einstein(var),'-o',label = 'Einstein theory')
ax.plot(var,debye(r/var),'-*',label='Debye theory')
ax.plot(var,[1]*len(var),'--',label = 'Dulong petit')
ax.set_title('Specific heat of a solid')
ax.set_xlabel('T/Te')
ax.set_ylim([0,1.1])
plt.legend(loc = 4)
plt.show()





