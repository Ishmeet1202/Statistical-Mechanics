import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

#FOR FERMI DIRAC 
# Here I am considering every constants in Electron volt

print("Reference - Garg Bansal Ghosh")
k=8.17e-5
T=[2000,12000]
m=0.511e-6
h=4.13e-15

def dn_de(T,E,meu,j):
    C=2*np.pi*(((2*m)**(3/2))/(h**3))
    if j==True:
        func=C * (E**(1/2))/(np.exp((E-meu)/(k*T))+1)   # Fermi- Dirac
    if j==False:
        func=C * (E**(1/2))/(np.exp((E-meu)/(k*T))-1)  # Bose- Eisntein
    return func

e=np.linspace(0,10,500)
print(e)

for i in T:    
    a=dn_de(i,e,5.49,True)
    plt.plot(e,a,label="T="+str(i))
    print(a)
plt.legend()
plt.xlabel("Energy")
plt.ylabel("dN/dE")
plt.grid()
plt.title("Fermi- Dirac(Non Relativistic)")
plt.show()




a=dn_de(2000,e,-5.49,False)
plt.plot(e,a,label="T=Low")
plt.legend()
plt.xlabel("Energy")
plt.ylabel("dN/dE")
plt.grid()
plt.title("Bose - Einstein (Non Relativistic)")
plt.show()

a=dn_de(10000,e,-5.49,False)
plt.plot(e,a,label="T=High")
plt.legend()
plt.xlabel("Energy")
plt.ylabel("dN/dE")
plt.grid()
plt.title("Bose - Einstein (Non Relativistic)")
plt.show()


def N(T,meu):
    func=lambda E: dn_de(T,E,meu,True)
    b=quad(func,0,1)
    return b

Y=[]
T=np.linspace(1000,5000,500)
for i in  T:
    Y.append(N(i,1)[0])
plt.scatter(T,Y )
plt.title("N vs T(Fermi- Dirac)")
plt.xlabel("Temperature")
plt.ylabel("N")
plt.grid()
plt.show()

# def N(T,meu):
#     func=lambda E: dn_de(T,E,meu,False)
#     b=quad(func,0,1)
#     return b

# Y=[]
# T=np.linspace(1000,1e5,500)
# for i in  T:
#     Y.append(N(i,-1)[0])
# plt.scatter(T,Y )
# plt.title("N vs T(Bose-Einstein)")
# plt.xlabel("Temperature")
# plt.ylabel("N")
# plt.grid()
# plt.show()

#For relativistic case

def dn_de_1(T,E,meu,j):
    C=2*np.pi*(((2*m)**(3/2))/(h**3))
    if j==True:
        func=C * (E**(2))/(np.exp((E-meu)/(k*T))+1)   # Fermi- Dirac
    if j==False:
        func=C * (E**(2))/(np.exp((E-meu)/(k*T))-1)  # Bose- Eisntein
    return func

e=np.linspace(0,10,500)
T=[2000,10000]
for i in T:    
    a=dn_de_1(i,e,5.49,True)
    plt.plot(e,a,label="T="+str(i))
plt.legend()
plt.xlabel("Energy")
plt.ylabel("dN/dE")
plt.grid()
plt.title("Fermi- Dirac (Relativistic)")
plt.show()

a=dn_de_1(2000,e,-5.49,False)
plt.plot(e,a,label="T=Low")
plt.legend()
plt.xlabel("Energy")
plt.ylabel("dN/dE")
plt.grid()
plt.title("Bose- Einstein (Relativistic)")
plt.show()

a=dn_de_1(10000,e,-5.49,False)
plt.plot(e,a,label="T=High")
plt.legend()
plt.xlabel("Energy")
plt.ylabel("dN/dE")
plt.grid()
plt.title("Bose- Einstein (Relativistic)")
plt.show()


# Energy Function
def Energy(T,E,meu,j):
    C=2*np.pi*(((2*m)**(3/2))/(h**3))
    if j==True:
        func=C * (E**(3/2))/(np.exp((E-meu)/(k*T))+1)   # Fermi- Dirac
    if j==False:
        func=C * (E**(3/2))/(np.exp((E-meu)/(k*T))-1)  # Bose- Eisntein
    return func
def U(T,meu):
    func=lambda E: Energy(T,E,meu,True)
    b=quad(func,0,10)
    return b

Y=[]
T=np.linspace(2000,1e6,500)
for i in  T:
    Y.append(U(i,5.49)[0])
plt.plot(T,Y )
plt.title("N.R. Energy Density(Fermi-Dirac)")
plt.xlabel("Temperature")
plt.ylabel("E/V")
plt.grid()
plt.show()

def U(T,meu):
    func=lambda E: Energy(T,E,meu,False)
    b=quad(func,0,10)
    return b

Y=[]
T=np.linspace(2000,1e7,500)
for i in  T:
    Y.append(U(i,-5.49)[0])
plt.plot(T,Y )
plt.title("N.R. Energy Density(Bose Einstein)")
plt.xlabel("Temperature")
plt.ylabel("E/V")
plt.grid()
plt.show()

def Energy(T,E,meu,j):
    C=2*np.pi*(((2*m)**(3/2))/(h**3))
    if j==True:
        func=C * (E**(3))/(np.exp((E-meu)/(k*T))+1)   # Fermi- Dirac
    if j==False:
        func=C * (E**(3))/(np.exp((E-meu)/(k*T))-1)  # Bose- Eisntein
    return func

def U(T,meu):
    func=lambda E: Energy(T,E,meu,True)
    b=quad(func,0,10)
    return b

Y=[]
T=np.linspace(2000,1e7,500)
for i in  T:
    Y.append(U(i,5.49)[0])
plt.plot(T,Y )
plt.title("Relativistic Energy Density(Fermi Dirac))")
plt.xlabel("Temperature")
plt.ylabel("E/V")
plt.grid()
plt.show()


Y=[]
T=np.linspace(2000,1e7,500)
for i in  T:
    Y.append(U(i,-5.49)[0])
plt.plot(T,Y )
plt.title("Relativistic Energy Density(Bose- Einstien)")
plt.xlabel("Temperature")
plt.ylabel("E/V")
plt.grid()
plt.show()