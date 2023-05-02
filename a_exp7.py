import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.stats import linregress
from numpy import diff

akarsh=2e-18
k=1.38e-23
N=6.02e23

def part_function(v,t):

    def function(n,v,t):
        return (n**2)*np.exp(-((n**2)*akarsh)/ (t*v**(2/3)))
    
    return (np.pi/2)*quad(function,1,1e11,(v,t))[0]


v=np.linspace(2e-2,5e-2,50)
t=np.linspace(150,450,50)

def matrix(v,t):
    result=np.zeros(((len(v)),(len(t))))
    result_log=np.zeros(((len(v)),(len(t))))
    for i in range(len(v)):
        for j in range(len(t)):
            result[i][j]=(part_function(v[i],t[j]))
            result_log[i][j]=np.log(part_function(v[i],t[j]))
    return result,result_log


slope=linregress( np.log(v),matrix(v,t)[1][:,0])
print("slope with v on x axis",slope[0])
for i in range(3):
    plt.plot(np.log(v),matrix(v,t)[1][:,i],label="t=" +str(t[i]))
plt.grid()
plt.legend()
plt.xlabel("Log volume")
plt.ylabel("Log Z")
plt.title("Z vs V")
plt.show()

slope=linregress( np.log(t),matrix(v,t)[1][0,:])
print("slope with t on x axis",slope[0])

for i in range(3):
    plt.plot(np.log(t), matrix(v,t)[1][i,:] , label="v=" +str(v[i]))
plt.legend()
plt.grid()
plt.xlabel("Log Temperature")
plt.ylabel("Log Z")
plt.title("Z vs T")
plt.show()


derivative= (diff(matrix(v,t)[1][:,0]))/(diff(v))

def pressure(t,d_t):
    return N*k*t *d_t
t1=[150,250,350,450]
for i in t1:
    plt.plot(v[:-1], pressure(i,derivative) , label='T=' +str(i))
plt.grid()
plt.legend()
plt.xlabel("Volume")
plt.ylabel("Pressure")
plt.title("P Vs V")
plt.show()

derivative_1=(diff(matrix(v,t)[1][0,:]))/(diff(t))

U=[]
def energy_density(temp , d_t):
    U.append(k*(np.array((temp))**2)*(d_t))
    return k*((np.array(temp))**2)*(d_t) ,U 

plt.plot(t[:-1], energy_density(t[:-1],derivative_1)[0])
plt.grid()
plt.xlabel("Temperature")
plt.ylabel("Energy ")
plt.title("Energy Vs T")
plt.show()

slope=linregress( t[:-1], energy_density(t[:-1],derivative_1)[0])
C_v =slope[0] 
print("C_v =" , C_v)

def varriance(t):
    return k*(t)**2* C_v

plt.plot(t, varriance(t))
plt.grid()
plt.xlabel("Temperature")
plt.ylabel(" Varriance")
plt.title("varriance Vs T")
plt.show()

def entropy (u, t):
    return (u*N)/(t) + (N*k*((matrix(v,t))[1][0,:] - np.log(N) +1))

plt.plot(t[:-1], entropy(np.array(U),t[:-1])[0])
plt.grid()
plt.xlabel("Temperature")
plt.ylabel("S")
plt.title("S Vs T")
plt.show()





    