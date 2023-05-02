import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import stats

k = 8.61e-5
N = 1

def partition(lev,T):
    e = np.arange(0,lev,1)
    g = np.array([1]*lev)
    fun = np.exp(-e/(k*T))
    
    z = g.dot(fun) 
    ni_n = (g*fun)/z
    u = (np.array(ni_n)*N).dot(e)
    s = (N*k*np.log(z/N)) + (u/T) + (N*k)
    f = -N*k*T*np.log(z)

    return z,ni_n[0],ni_n[1],ni_n[2],u,s,f

partition = np.vectorize(partition,otypes = [np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray])

T1 = np.linspace(1e-5,5e3,100)
T2 =  np.linspace(5e3,1e5,100)

z1,ni11,ni12,ni13,u1,s1,f1 = partition(3,T1)
z2,ni21,ni22,ni23,u2,s2,f2 = partition(3,T2)

def plot(x1,x2,y1,y2,title,ylabel,key = None,y3 = None,y4=None,y5 = None,y6 = None):
    fig,ax = plt.subplots(1,2)
    plt.suptitle(title)
    ax[0].grid()
    ax[1].grid()
    ax[0].plot(x1,y1,'-o',label ='Low temperature')
    ax[1].plot(x2,y2,'-v',label = 'High temperature')
    if key == 1:
       ax[0].plot(x1,y3,'r-o',label ='Low temperature')
       ax[1].plot(x2,y4,'r-v',label = 'High temperature')
       ax[0].plot(x1,y5,'g-o',label ='Low temperature')
       ax[1].plot(x2,y6,'g-v',label = 'High temperature')
       ax[1].plot(x2,[1/3]*len(x2))
    ax[0].set_xlabel('T')
    ax[1].set_xlabel('T') 
    ax[0].set_ylabel(ylabel)
    
    plt.show()


res = stats.linregress(T2,f2)

print('entropy', -res[0])


#io = np.where(abs(s1+res[0]) < (-res)/200)

#print(T1[io])

plot(T1,T2,z1,z2,'Partition function','z')
plot(T1,T2,ni11,ni21,'fractionl particles','ni/n',1,ni12,ni22,ni13,ni23)

plot(T1,T2,u1,u2,'Total energy','U')
plot(T1,T2,s1,s2,'Total entropy','S')
plot(T1,T2,f1,f2,'Helmholtz free energy','F')



      
    
    
    
     
    

