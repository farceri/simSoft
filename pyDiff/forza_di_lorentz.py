import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.stats import norm
import sys
import sympy as sp



dt=0.001
N= 100
gamma=0.05 #attenzione a sceglierla perchè ti incasina tutto
temperature=1.2
kB=1.0
mass=1.0 #dovrebbe essere 1.0, ma ho problemi di overflow
sigma=np.sqrt(2*mass*gamma*kB*temperature)

v_0x=13
v_0y=20
v_0=np.sqrt((v_0x)**2+(v_0y)**2)
x_0=0

#VERIFICA PDF RISPETTO A D_VV (GAUSSIANA)
x_t=np.linspace(0,2000,1000)
t_fiss=450
def x_mediat(x_0,v_0,gamma,t_fiss):
        return( x_0+(v_0/gamma)*(1-np.exp(-gamma*t_fiss))
        )


def MSDt(v_0,sigma,gamma,mass,t_fiss):
        return(
            ((sigma/gamma*t_fiss)**2)*t_fiss + (((sigma)**2)/(2*(mass**2)*(gamma**3)))*(4*np.exp(-gamma*t_fiss)-np.exp(-2*gamma*t_fiss)-3) + (v_0**2)*((1/gamma)*(np.exp(-gamma*t_fiss)-1))**2
        )

D_vv=((sigma)**2)/2

varianza=4*D_vv*t_fiss
if t_fiss >=200:
       sigma1=varianza
else:
       sigma1=MSDt(v_0,sigma,gamma,mass,t_fiss)

def pdf(sigma1,x_t,v_0,x_0,gamma,t_fiss):
       media = x_mediat(x_0,v_0,gamma,t_fiss)
       return ((1/np.sqrt(2*np.pi*(sigma1**2)))*np.exp(-((x_t-media)**2)/(2*sigma1**2)))


print("Media = ", x_mediat(x_0,v_0,gamma,t_fiss))
print("Sigma = ", sigma1)

dati = np.random.normal(loc=x_mediat(x_0,v_0,gamma,t_fiss), scale=sigma1, size=10000)
conteggi, bin_edges = np.histogram(dati, bins=50, density=True)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
plt.plot(bin_centers, conteggi, label="Istogramma normalizzato")
x = x_t
plt.plot(x, norm.pdf(x, loc=x_mediat(x_0,v_0,gamma,t_fiss), scale=sigma1), label="PDF teorica", linestyle='dashed')
plt.xlabel("r")
plt.ylabel("Gaussiana")
plt.legend()
plt.grid(False)
plt.show()

#INTEGRALE DELL'ENERGIA 
N = 1000        
dt = 0.001        
time = np.arange(0, N*dt, dt)
v0=np.sqrt(kB *temperature /mass)
v_0=np.array([v0,v0])
r_0=np.array([0,0])
r=np.empty(shape=(N,2))
v=np.empty(shape=(N,2))
E = np.zeros(N)  
msd=np.zeros(N)  
r[0] = r_0     
v[0] = v_0

#integrazione verlet
for i in range(1, N):
       noise=np.sqrt(2 * kB *temperature *gamma /dt)*np.random.normal(0,1/dt)
       forza_trascinamento=-gamma * v[i-1]
       langevin_force=forza_trascinamento+noise
       kinetic_energy = 0.5 *mass *((np.linalg.norm(v[i-1]))** 2)
       compute_msd=np.mean((r[i-1])**2)
       v[i] = v[i-1] + 0.5*langevin_force/ mass * dt
       r[i] = r[i-1] + v[i-1]*dt 
       E[i] =kinetic_energy
       msd[i]=compute_msd
plt.figure(figsize=(10, 6))
# Grafico della posizione rispetto al tempo
plt.subplot(2, 1, 1)
plt.plot(time, msd)
plt.xlabel('t')
plt.ylabel('msd')
# Grafico dell'energia rispetto al tempo
plt.subplot(2, 1, 2)
plt.plot(time, E)
plt.xlabel('t')
plt.ylabel('k')
plt.title('k')
plt.tight_layout()
plt.show()

'''
#INTEGRALE DELL'ENERGIA 
N = 1000        
dt = 0.01        
time = np.arange(0, N*dt, dt)
v0=np.sqrt(kB *temperature /mass)
v_0=np.array([v0,v0])
r_0=np.array([0,0])
r=np.empty(shape=(N,2))
v=np.empty(shape=(N,2))
E = np.zeros(N)  
msd=np.zeros(N)  
r[0] = r_0     
v[0] = v_0
direzione=np.array([0,0])


#integrazione verlet
for i in range(1, N):
       asse = random.randint(0, len(direzione) - 1)
       direzione[asse] += 1
       noise=np.sqrt(2 * kB *temperature *gamma /time[i-1])*direzione[asse]
       forza_stokes=-gamma * v[i-1]
       langevin_force=forza_stokes+noise
       kinetic_energy = 0.5 *mass *((np.linalg.norm(v[i-1]))** 2)
       compute_msd=np.mean((r[i-1])**2)
       r[i] = r[i-1] + v[i-1]*dt 
       v[i] = v[i-1] + 0.5*langevin_force/ mass * dt
       E[i] =kinetic_energy
       msd[i]=compute_msd
plt.figure(figsize=(10, 6))
# Grafico della posizione rispetto al tempo
plt.subplot(2, 1, 1)
plt.plot(time, msd)
plt.xlabel('t')
plt.ylabel('msd')
# Grafico dell'energia rispetto al tempo
plt.subplot(2, 1, 2)
plt.plot(time, E)
plt.xlabel('t')
plt.ylabel('k')
plt.title('k')
plt.tight_layout()
plt.show()
'''


#FLUSSO


#FORZA DI LORENTZ
q=1
v_fl=np.array([1,2,0])
B=np.array([0,0,20])
forza_di_lorentz=np.cross(q*v_fl,B)
print('La forza di LORENTZ è:',forza_di_lorentz)

#INTEGRALE ENERGIA CON FORZA DI LORENTZ 
N = 1000        
dt = 0.001    
time = np.arange(0, N*dt, dt)
v0= np.sqrt(kB *temperature /mass)
v_0=np.array([v0,v0,v0])
r_0=np.array([0,0,0])
r=np.empty(shape=(N,3))
v=np.empty(shape=(N,3))
E=np.zeros(N)
msd=np.zeros(N)
v[0]=v_0
r[0]=r_0
direzione=np.array([0,0,0])

for i in range (1,N):
        noise=np.sqrt(2 * kB *temperature *gamma /dt)*np.random.normal(0,1/dt)
        forza_trascinamento=-gamma *v[i-1]
        forza_di_lorentz=np.cross(q*v[i-1],B)
        forza_langevin=forza_trascinamento+noise+forza_di_lorentz
        kinetic_energy = 0.5 *mass *((np.linalg.norm(v[i-1]))** 2)
        compute_msd=np.mean(np.sum(np.linalg.norm(r[i-1]))**2)
        r[i] = r[i-1] + v[i-1]*dt 
        v[i] = v[i-1] + 0.5*forza_langevin/ mass * dt
        E[i] =kinetic_energy
        msd[i]=compute_msd

plt.figure(figsize=(10, 6))
# Grafico della posizione rispetto al tempo
plt.subplot(2, 1, 1)
plt.plot(time, msd)
plt.xlabel('t')
plt.ylabel('msd')
# Grafico dell'energia rispetto al tempo
plt.subplot(2, 1, 2)
plt.plot(time, E)
plt.xlabel('t')
plt.ylabel('E')
plt.tight_layout()
plt.show()


'''
#FORZA DI LORENTZ
q=1
v_fl=np.array([1,2,0])
B=np.array([0,0,5])
forza_di_lorentz=np.cross(q*v_fl,B)
print('La forza di LORENTZ è:',forza_di_lorentz)

#FORZA DI LANGEVIN
noise=-gamma * v_fl+np.sqrt(2 * kB *temperature *gamma /dt)
forza_langevin=noise+forza_di_lorentz
print('La forza di Langevin è:', forza_langevin)

#INTEGRALE ENERGIA CON FORZA DI LORENTZ 2
N = 1000        
dt = 15     
time = np.arange(0, N*dt, dt)
v0= np.sqrt(kB *temperature /mass)
v_0=np.array([v0,v0,v0])
r_0=np.array([0,0,0])
r=np.empty(shape=(N,3))
v=np.empty(shape=(N,3))
E=np.zeros(N)
msd=np.zeros(N)
v[0]=v_0
r[0]=r_0
direzione=np.array([0,0,0])

for i in range (1,N):
        asse = random.randint(0, len(direzione) - 2)
        direzione[asse] += 1
        noise=np.sqrt(2 * kB *temperature *gamma /time[i-1])*direzione[asse]
        forza_stokes=-gamma *v[i-1]
        forza_di_lorentz=np.cross(q*v[i-1],B)
        forza_langevin=forza_stokes+noise+forza_di_lorentz
        kinetic_energy = 0.5 *mass *((np.linalg.norm(v[i-1]))** 2)
        compute_msd=np.mean(np.sum(np.linalg.norm(r[i-1]))**2)
        r[i] = r[i-1] + v[i-1]*dt 
        v[i] = v[i-1] + 0.5*forza_langevin/ mass * dt
        E[i] =kinetic_energy
        msd[i]=compute_msd

plt.figure(figsize=(10, 6))
# Grafico della posizione rispetto al tempo
plt.subplot(2, 1, 1)
plt.plot(time, msd)
plt.xlabel('t')
plt.ylabel('msd')
# Grafico dell'energia rispetto al tempo
plt.subplot(2, 1, 2)
plt.plot(time, E)
plt.xlabel('t')
plt.ylabel('E')
plt.tight_layout()
plt.show()
'''

#COEFFICIENTE DI DIFFUSIONE
k=(q*B[2])/gamma
D_0=(kB*temperature)/(mass*gamma)
m=np.array([[1,k],[-k,1]])
D=(D_0/(1+k**2))*m
print('D=',D)




#TEMPO DI DIFFUSIONE
n=100
t=np.arange(1,100,2)
L=5 #sistem size
diffusivity=D_0/(1+(k**2))
print(diffusivity)
time_of_diffusion=(D_0*t)/((1+(k**2))*(L**2))
xval=t
yval=time_of_diffusion
plt.xlabel('t',fontsize=15)
plt.ylabel('Tempo di diffusione',fontsize=15)
plt.plot(xval,yval)
plt.show()

'''
#D per tempo di diffusione al variare del campo magnetico
B=np.array([0,0,5])
num_step = 4
for step in range(num_step):
    B[2] += 1
k=(q*B[2])/gamma
D_0=(kB*temperature)/(mass*gamma)
t=1000
L=1
time_of_diffusion=(D_0*t)/((1+(k**2))*(L**2))
xval=k
yval=time_of_diffusion
plt.xlabel('k',fontsize=15)
plt.ylabel('Tempo di diffusione',fontsize=15)
plt.plot(xval,yval)
plt.show()
'''
    

#ESEMPIO PER CALCOLARE IL GRADIENTE DI UN VETTORE
f=np.array([1,2,4,7,11,16])
grad_f1=np.gradient(f,1)
print(grad_f1)
grad_f2=np.gradient(f,2)
print(grad_f2)

#gradiente coefficiente di diffusione
grad_D=np.gradient(D,1)
print(grad_D)

'''
#FLUSSO
L=1
E=5
x_t=sp.symbols('x(t)')
prob=sp.pdf(sigma1,x_t,v_0,x_0,gamma,t_fiss)
j_fl=D*((q*E)/(kB*temperature))
print(j_fl)

'''




