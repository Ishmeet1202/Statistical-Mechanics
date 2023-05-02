import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as int1
from scipy.stats import linregress

def F_p(x):
    return (x**3)/(np.exp(x)-1)

x=np.linspace(1e-8,12,151)
X=np.array(x)
P=np.array(F_p(x))
m =  np.max(F_p(x))
t = F_p(x)
i_val = np.where(t == m)

plt.scatter(x,F_p(x), marker='*')
plt.grid()
plt.xlabel("x")
plt.ylabel("F_p")
plt.title("F_p")
plt.show()
Area=int1.quad(F_p,0,12)[0]
print("Total area:-" , Area)
print("X Median:- " , Area/2)

h=6.62e-34
c=3e8
k=1.38*10**(-23)

print(" Value of Weins Displacement Constant for x_m:-",h*c/(k*Area/2))

# I = int1.quad(F_p ,0, 40)[0]
I_1=int1.simps(P,X)
print("The value of I_p in Stefan- Boltzman",I_1)
# print("The value of I_p in Stefan- Boltzman" , I)

def U(T):
    return (np.pi**4)/(15) * 8*(np.pi)*((k*T)**4/(h*c)**3)

T=np.arange(100,10000,500)
values=[]
for i in T:
    value=U(i)*(15/(np.pi)**4)
    values.append(value)
# print(values )

def radiant_flux(T,C):
    return (c/4)*((np.pi)**4/(15)) *C

for i in T:
    value=U(i)*(15/(np.pi)**4)
    plt.scatter(i,radiant_flux(i,value) , label='T =' + str(i))
plt.legend()
plt.xlabel("Temperature")
plt.ylabel("F(T)")
plt.grid()
plt.show()

temp=[]
for i in T:
    value=U(i)*(15/(np.pi)**4)
    plt.scatter(np.log(i),np.log(radiant_flux(i,value)))
    temp.append(radiant_flux(i,value))
slope=linregress(np.log(T) , np.log(temp))
print("slope",slope[0])
print("intercept",(slope[1]))
plt.legend()
plt.xlabel("Temperature")
plt.ylabel("F(T)")
plt.title("Log plot of F(T) and Temperature")
plt.grid()
plt.show()
print ("Stefan- Boltzman constant value :-" , np.exp(slope[1]))
