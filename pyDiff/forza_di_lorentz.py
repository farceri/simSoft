import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.stats import norm
import sys
import sympy as sp
from scipy.signal import butter, lfilter


dt=0.001
N= 100
gamma=0.05 
temperature=1.2
kB=1.0
mass=1.0 
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


#integrale energia senza rumore
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
       forza_trascinamento=-gamma * v[i-1]
       kinetic_energy = 0.5 *mass *((np.linalg.norm(v[i-1]))** 2)
       compute_msd=np.mean((r[i-1])**2)
       v[i] = v[i-1] + 0.5*forza_trascinamento/ mass * dt
       r[i] = r[i-1] + v[i-1]*dt
       E[i] =kinetic_energy
       msd[i]=compute_msd
       

plt.figure(figsize=(10, 6))
# Grafico delle componenti della posizione rispetto al tempo
plt.subplot(2, 1, 1)
plt.plot(time, r[:,0], label='componente x')
plt.plot(time, r[:,1], label="componente y")
plt.xlabel('t')
plt.legend()
plt.ylabel('posizione')
# Grafico delle componenti della velocità rispetto al tempo
plt.subplot(2, 1, 2)
plt.plot(time, v[:,0], label='componente x')
plt.plot(time, v[:,1], label="componente y")
plt.xlabel('t')
plt.ylabel('velocità')
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
# Grafico msd rispetto al tempo
plt.subplot(2, 1, 1)
plt.plot(time, msd)
plt.xlabel('t')
plt.ylabel('msd')
# Grafico dell'energia rispetto al tempo
plt.subplot(2, 1, 2)
plt.plot(time, E)
plt.xlabel('t')
plt.ylabel('k')
plt.tight_layout()
plt.show()



#RUMORE
N = 1000        
dt = 0.001        
time = np.arange(0, N*dt, dt)
rumore=np.zeros(shape=(N,2))
rumore[0]=0
for i in range(1, N):
       rumore_bianco=np.sqrt(2 * kB *temperature *gamma/dt)*np.array([np.random.normal(0,1),np.random.normal(0,1)])
       rumore[i]=rumore_bianco*dt
plt.figure(figsize=(10, 6))
plt.plot(time, rumore[:,0], label='componente x')
plt.plot(time, rumore[:,1], label="componente y")
plt.xlabel('t')
plt.ylabel('noise')
plt.show()



#INTEGRALE DELL'ENERGIA con il rumore(bidimensionale)
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
       noise=np.sqrt(2 * kB *temperature *gamma /dt)*np.array([np.random.normal(0,1),np.random.normal(0,1)])
       forza_trascinamento=-gamma * v[i-1]
       langevin_force=forza_trascinamento+noise
       kinetic_energy = 0.5 *mass *((np.linalg.norm(v[i-1]))** 2)
       compute_msd=np.mean((r[i-1])**2)
       v[i] = v[i-1] + 0.5*langevin_force/ mass * dt
       r[i] = r[i-1] + v[i-1]*dt +noise*dt
       E[i] =kinetic_energy
       msd[i]=compute_msd
       

plt.figure(figsize=(10, 6))
# Grafico delle componenti della posizione rispetto al tempo
plt.subplot(2, 1, 1)
plt.plot(time, r[:,0], label='componente x')
plt.plot(time, r[:,1], label="componente y")
plt.xlabel('t')
plt.legend()
plt.ylabel('posizione')
# Grafico delle componenti della velocità rispetto al tempo
plt.subplot(2, 1, 2)
plt.plot(time, v[:,0], label='componente x')
plt.plot(time, v[:,1], label="componente y")
plt.xlabel('t')
plt.ylabel('velocità')
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
# Grafico msd rispetto al tempo
plt.subplot(2, 1, 1)
plt.plot(time, msd)
plt.xlabel('t')
plt.ylabel('msd')
# Grafico dell'energia rispetto al tempo
plt.subplot(2, 1, 2)
plt.plot(time, E)
plt.xlabel('t')
plt.ylabel('k')
plt.tight_layout()
plt.show()


#integrale energia confronto con e senza rumore per componenti
N = 1000        
dt = 0.001        
time = np.arange(0, N*dt, dt)
v0=np.sqrt(kB *temperature /mass)
v_0=np.array([v0,v0])
r_0=np.array([0,0])
r=np.empty(shape=(N,2))
r_s=np.empty(shape=(N,2))
v=np.empty(shape=(N,2))
v_s=np.empty(shape=(N,2))
E = np.zeros(N)
E_s = np.zeros(N)  
msd=np.zeros(N)  
msd_s=np.zeros(N)  
r[0] = r_0  
r_s[0]=r_0   
v[0] = v_0
v_s[0]=v_0


#integrazione verlet
for i in range(1, N):
       noise=np.sqrt(2 * kB *temperature *gamma /dt)*np.array([np.random.normal(0,1),np.random.normal(0,1)])
       forza_trascinamento=-gamma * v[i-1]
       langevin_force=forza_trascinamento+noise
       kinetic_energy = 0.5 *mass *((np.linalg.norm(v[i-1]))** 2)
       k_s=0.5 *mass *((np.linalg.norm(v_s[i-1]))** 2)
       compute_msd=np.mean((r[i-1])**2)
       c_msd_s=np.mean((r_s[i-1])**2)
       v[i] = v[i-1] + 0.5*langevin_force/ mass * dt
       v_s[i]=v_s[i-1] + 0.5*forza_trascinamento/ mass * dt
       r[i] = r[i-1] + v[i-1]*dt +noise*dt
       r_s[i] = r_s[i-1] + v_s[i-1]*dt 
       E[i] =kinetic_energy
       E_s[i] =k_s
       msd[i]=compute_msd
       msd_s[i]=c_msd_s

plt.figure(figsize=(10, 6))
# Grafico delle componenti della posizione componente x
plt.subplot(2, 1, 1)
plt.plot(time, r_s[:,0], label='no noise')
plt.plot(time, r[:,0], label="with noise")
plt.xlabel('t')
plt.legend()
plt.ylabel('x')
# Grafico delle componenti della posizione componente y
plt.subplot(2, 1, 2)
plt.plot(time, r_s[:,1], label='no noise')
plt.plot(time, r[:,1], label="with noise")
plt.xlabel('t')
plt.legend()
plt.ylabel('y')
plt.subplots_adjust(hspace=0.3) 
plt.show()


plt.figure(figsize=(10, 6))
# Grafico delle componenti della velocità componente x
plt.subplot(2, 1, 1)
plt.plot(time, v_s[:,0], label='no noise')
plt.plot(time, v[:,0], label="with noise")
plt.xlabel('t')
plt.legend()
plt.ylabel('v_x')
# Grafico delle componenti della velocità componente y
plt.subplot(2, 1, 2)
plt.plot(time, v_s[:,1], label='no noise')
plt.plot(time, v[:,1], label="with noise")
plt.xlabel('t')
plt.legend()
plt.ylabel('v_y')
plt.subplots_adjust(hspace=0.3)
plt.show()


# Grafico msd rispetto al tempo
plt.figure(figsize=(8, 6))
plt.plot(time, msd_s, label='no noise')
plt.plot(time,msd,label="with noise")
plt.xlabel('t')
plt.ylabel('msd')
plt.legend()
plt.show()


# Grafico dell'energia rispetto al tempo
plt.figure(figsize=(8, 6))
plt.plot(time, E_s, label='no noise')
plt.plot(time,E,label="with noise")
plt.xlabel('t')
plt.ylabel('E')
plt.legend()
plt.show()

#lg eq m=0
N = 1000        
dt = 0.001
gamma_1=200       
time = np.arange(0, N*dt, dt)
v0=np.sqrt(kB *temperature /mass)
r_0=np.array([0,0])
r=np.empty(shape=(N,2))
msd=np.zeros(N)  
r[0] = r_0     

#integrazione verlet
for i in range(1, N):
       noise=np.sqrt(2 * kB *temperature *gamma_1 /dt)*np.array([np.random.normal(0,1),np.random.normal(0,1)])
       compute_msd=np.mean((r[i-1])**2)
       r[i] = r[i-1] +(noise*dt)/gamma_1
       msd[i]=compute_msd

plt.figure(figsize=(10, 6))
# Grafico delle componenti della posizione rispetto al tempo
plt.subplot(2, 1, 1)
plt.plot(time, r[:,0], label='componente x')
plt.plot(time, r[:,1], label="componente y")
plt.xlabel('t')
plt.legend()
plt.ylabel('posizione')

plt.subplot(2, 1, 2)
plt.plot(time, msd)
plt.xlabel('t')
plt.ylabel('msd')
plt.show()

#FORZA DI LORENTZ
q=1
v_fl=np.array([1,2,0])
B=np.array([0,0,20])
forza_di_lorentz=np.cross(q*v_fl,B)
print('La forza di LORENTZ è:',forza_di_lorentz)

#integrale energia con forza di Lorentz
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
       noise=np.sqrt(2 * kB *temperature *gamma)*np.array([[np.random.normal(0,1),np.random.normal(0,1)]])
       G=np.array([[gamma,-q*B[2]],[q*B[2],gamma]])
       forza_trascinamento=np.dot(-G , v[i-1])
       kinetic_energy = 0.5 *mass *((np.linalg.norm(v[i-1]))** 2)
       compute_msd=np.mean((r[i-1])**2)
       v[i] = v[i-1] + 0.5*forza_trascinamento/ mass * dt + noise
       r[i] = r[i-1] + v[i-1]*dt +noise
       E[i] =kinetic_energy
       msd[i]=compute_msd
       

plt.figure(figsize=(10, 6))
# Grafico delle componenti della posizione rispetto al tempo
plt.subplot(2, 1, 1)
plt.plot(time, r[:,0], label='componente x')
plt.plot(time, r[:,1], label="componente y")
plt.xlabel('t')
plt.legend()
plt.ylabel('posizione')
# Grafico delle componenti della velocità rispetto al tempo
plt.subplot(2, 1, 2)
plt.plot(time, v[:,0], label='componente x')
plt.plot(time, v[:,1], label="componente y")
plt.xlabel('t')
plt.ylabel('velocità')
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
# Grafico msd rispetto al tempo
plt.subplot(2, 1, 1)
plt.plot(time, msd)
plt.xlabel('t')
plt.ylabel('msd')
# Grafico dell'energia rispetto al tempo
plt.subplot(2, 1, 2)
plt.plot(time, E)
plt.xlabel('t')
plt.ylabel('k')
plt.tight_layout()
plt.show()

#integrale energia confronto con e senza forza di lorentz per componenti
N = 1000        
dt = 0.001        
time = np.arange(0, N*dt, dt)
v0=np.sqrt(kB *temperature /mass)
v_0=np.array([v0,v0])
r_0=np.array([0,0])
r=np.empty(shape=(N,2))
v=np.empty(shape=(N,2))
r_fl=np.empty(shape=(N,2))
v_fl=np.empty(shape=(N,2))
r=np.empty(shape=(N,2))
v=np.empty(shape=(N,2))
E_fl = np.zeros(N)  
msd_fl=np.zeros(N)  
E = np.zeros(N)  
msd=np.zeros(N)  
r[0] = r_0     
v[0] = v_0
r_fl[0] = r_0     
v_fl[0] = v_0


#integrazione verlet
for i in range(1, N):
       noise=np.sqrt(2 * kB *temperature *gamma)*np.array([[np.random.normal(0,1),np.random.normal(0,1)]])
       f_t=-gamma * v[i-1]
       G=np.array([[gamma,-q*B[2]],[q*B[2],gamma]])
       forza_trascinamento=np.dot(-G , v_fl[i-1])
       kinetic_energy = 0.5 *mass *((np.linalg.norm(v_fl[i-1]))** 2)
       k=0.5 *mass *((np.linalg.norm(v[i-1]))** 2)
       compute_msd=np.mean((r_fl[i-1])**2)
       c_msd=np.mean((r[i-1])**2)
       v[i] = v[i-1] + 0.5*f_t/ mass * dt + noise
       r[i] = r[i-1] + v[i-1]*dt +noise
       v_fl[i] = v_fl[i-1] + 0.5*forza_trascinamento/ mass * dt + noise
       r_fl[i] = r_fl[i-1] + v_fl[i-1]*dt +noise
       E_fl[i] =kinetic_energy
       msd_fl[i]=compute_msd
       E[i] =k
       msd[i]=c_msd

plt.figure(figsize=(10, 6))
# Grafico delle componenti della posizione componente x
plt.subplot(2, 1, 1)
plt.plot(time, r[:,0], label='no lorentz_force')
plt.plot(time, r_fl[:,0], label="with lorentz_force")
plt.xlabel('t')
plt.legend()
plt.ylabel('x')
# Grafico delle componenti della posizione componente y
plt.subplot(2, 1, 2)
plt.plot(time, r[:,1], label='no lorentz_force')
plt.plot(time, r_fl[:,1], label="with lorentz_force")
plt.xlabel('t')
plt.legend()
plt.ylabel('y')
plt.subplots_adjust(hspace=0.3) 
plt.show()

# Grafico msd rispetto al tempo
plt.figure(figsize=(8, 6))
plt.plot(time, msd_fl, label='no lorentz_force')
plt.plot(time,msd,label="with lorentz_force")
plt.xlabel('t')
plt.ylabel('msd')
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
# Grafico delle componenti della velocità componente x
plt.subplot(2, 1, 1)
plt.plot(time, v[:,0], label='no lorentz_force')
plt.plot(time, v_fl[:,0], label="with lorentz_force")
plt.xlabel('t')
plt.legend()
plt.ylabel('v_x')
# Grafico delle componenti della velocità componente y
plt.subplot(2, 1, 2)
plt.plot(time, v[:,1], label='no lorentz_force')
plt.plot(time, v_fl[:,1], label="with lorentz_force")
plt.xlabel('t')
plt.legend()
plt.ylabel('v_y')
plt.subplots_adjust(hspace=0.3)
plt.show()

# Grafico dell'energia rispetto al tempo
plt.figure(figsize=(8, 6))
plt.plot(time, E, label='no lorentz_force')
plt.plot(time,E_fl,label="with lorentz_force")
plt.xlabel('t')
plt.ylabel('E')
plt.legend()
plt.show()


#FP eq corrispondente
G_s=np.linalg.inv(G)*np.matrix.transpose(np.linalg.inv(G))*gamma


#CONFRONTO TRA TUTTI E TRE (LOGARITMICO)
N = 1000        
dt = 0.001        
time = np.arange(0, N*dt, dt)
v0=np.sqrt(kB *temperature /mass)
v_0=np.array([v0,v0])
r_0=np.array([0,0])
r=np.empty(shape=(N,2))
v=np.empty(shape=(N,2))
r_fl=np.empty(shape=(N,2))
v_fl=np.empty(shape=(N,2))
r_sn=np.empty(shape=(N,2))
v_sn=np.empty(shape=(N,2))
E_fl = np.zeros(N)  
msd_fl=np.zeros(N)  
E = np.zeros(N)  
msd=np.zeros(N)  
E_sn = np.zeros(N)  
msd_sn=np.zeros(N)  
r[0] = r_0     
v[0] = v_0
r_fl[0] = r_0     
v_fl[0] = v_0
r_sn[0] = r_0     
v_sn[0] = v_0

#integrazione verlet
for i in range(1, N):
       noise=np.sqrt(2 * kB *temperature *gamma)*np.array([[np.random.normal(0,1),np.random.normal(0,1)]])
       f_t=-gamma * v[i-1]
       G=np.array([[gamma,-q*B[2]],[q*B[2],gamma]])
       forza_trascinamento=np.dot(-G , v_fl[i-1])
       kinetic_energy = 0.5 *mass *((np.linalg.norm(v_fl[i-1]))** 2)
       k=0.5 *mass *((np.linalg.norm(v[i-1]))** 2)
       k_sn=0.5 *mass *((np.linalg.norm(v_sn[i-1]))** 2)
       compute_msd=np.mean((r_fl[i-1])**2)
       compute_msd_sn=np.mean((r_sn[i-1])**2)
       c_msd=np.mean((r[i-1])**2)
       v[i] = v[i-1] + 0.5*f_t/ mass * dt + noise
       r[i] = r[i-1] + v[i-1]*dt +noise
       v_sn[i] = v_sn[i-1] + 0.5*f_t/ mass * dt 
       r_sn[i] = r_sn[i-1] + v_sn[i-1]*dt 
       v_fl[i] = v_fl[i-1] + 0.5*forza_trascinamento/ mass * dt + noise
       r_fl[i] = r_fl[i-1] + v_fl[i-1]*dt +noise
       E_fl[i] =kinetic_energy
       msd_fl[i]=compute_msd
       E[i] =k
       msd[i]=c_msd
       E_sn[i] =k_sn
       msd_sn[i]=compute_msd_sn


plt.figure(figsize=(10, 6))
# Grafico delle componenti della posizione componente x
plt.subplot(2, 1, 1)
plt.plot(time, r_sn[:,0], label="no noise")
plt.plot(time, r[:,0], label='no lorentz_force')
plt.plot(time, r_fl[:,0], label="with lorentz_force")
plt.xlabel('t')
plt.legend()
plt.ylabel('x')
# Grafico delle componenti della posizione componente y
plt.subplot(2, 1, 2)
plt.plot(time, r_sn[:,1], label="no noise")
plt.plot(time, r[:,1], label='no lorentz_force')
plt.plot(time, r_fl[:,1], label="with lorentz_force")
plt.xlabel('t')
plt.legend()
plt.ylabel('y')
plt.subplots_adjust(hspace=0.3) 
plt.show()       

plt.figure(figsize=(10, 6))
# Grafico delle componenti della posizione componente x
plt.subplot(2, 1, 1)
plt.yscale('log')
plt.plot(time, r_sn[:,0], label="no noise")
plt.plot(time, r[:,0], label='no lorentz_force')
plt.plot(time, r_fl[:,0], label="with lorentz_force")
plt.xlabel('t')
plt.legend()
plt.ylabel('x')
# Grafico delle componenti della posizione componente y
plt.subplot(2, 1, 2)
plt.yscale('log')
plt.plot(time, r_sn[:,1], label="no noise")
plt.plot(time, r[:,1], label='no lorentz_force')
plt.plot(time, r_fl[:,1], label="with lorentz_force")
plt.xlabel('t')
plt.legend()
plt.ylabel('y')
plt.subplots_adjust(hspace=0.3) 
plt.show()

# Grafico msd rispetto al tempo
plt.figure(figsize=(8, 6))
plt.plot(time,msd_sn,label="no noise")
plt.plot(time, msd_fl, label='no lorentz_force')
plt.plot(time,msd,label="with lorentz_force")
plt.xlabel('t')
plt.ylabel('msd')
plt.legend()
plt.show()

plt.figure(figsize=(8, 6))
plt.yscale('log')
plt.plot(time,msd_sn,label="no noise")
plt.plot(time, msd_fl, label='no lorentz_force')
plt.plot(time,msd,label="with lorentz_force")
plt.xlabel('t')
plt.ylabel('msd')
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
# Grafico delle componenti della velocità componente x
plt.subplot(2, 1, 1)
plt.plot(time, v_sn[:,0], label='no noise')
plt.plot(time, v[:,0], label='no lorentz_force')
plt.plot(time, v_fl[:,0], label="with lorentz_force")
plt.xlabel('t')
plt.legend()
plt.ylabel('v_x')
# Grafico delle componenti della velocità componente y
plt.subplot(2, 1, 2)
plt.plot(time, v_sn[:,1], label='no noise')
plt.plot(time, v[:,1], label='no lorentz_force')
plt.plot(time, v_fl[:,1], label="with lorentz_force")
plt.xlabel('t')
plt.legend()
plt.ylabel('v_y')
plt.subplots_adjust(hspace=0.3)
plt.show()

plt.figure(figsize=(10, 6))
# Grafico delle componenti della velocità componente x
plt.subplot(2, 1, 1)
plt.yscale('log')
plt.plot(time, v_sn[:,0], label='no noise')
plt.plot(time, v[:,0], label='no lorentz_force')
plt.plot(time, v_fl[:,0], label="with lorentz_force")
plt.xlabel('t')
plt.legend()
plt.ylabel('v_x')
# Grafico delle componenti della velocità componente y
plt.subplot(2, 1, 2)
plt.yscale('log')
plt.plot(time, v_sn[:,1], label='no noise')
plt.plot(time, v[:,1], label='no lorentz_force')
plt.plot(time, v_fl[:,1], label="with lorentz_force")
plt.xlabel('t')
plt.legend()
plt.ylabel('v_y')
plt.subplots_adjust(hspace=0.3)
plt.show()

# Grafico dell'energia rispetto al tempo
plt.figure(figsize=(8, 6))
plt.plot(time, E_sn, label='no noise')
plt.plot(time, E, label='no lorentz_force')
plt.plot(time,E_fl,label="with lorentz_force")
plt.xlabel('t')
plt.ylabel('E')
plt.legend()
plt.show()

# Grafico dell'energia rispetto al tempo
plt.figure(figsize=(8, 6))
plt.yscale('log')
plt.plot(time, E_sn, label='no noise')
plt.plot(time, E, label='no lorentz_force')
plt.plot(time,E_fl,label="with lorentz_force")
plt.xlabel('t')
plt.ylabel('E')
plt.legend()
plt.show()


#lg eq m=0
N = 1000        
dt = 0.001
gamma_1=50       
time = np.arange(0, N*dt, dt)
v0=np.sqrt(kB *temperature /mass)
r_0=np.array([0,0])
r=np.empty(shape=(N,2))
rumore=np.empty(shape=(N,2))
msd=np.zeros(N)  
r[0] = r_0     

#integrazione verlet
for i in range(1, N):
       G=np.array([[gamma,-q*B[2]],[q*B[2],gamma]])
       G_inv=np.linalg.inv(G)
       G_inv_t=np.transpose(G_inv)
       eta=np.array([np.random.normal(0,1),np.random.normal(0,1)]) 
       noise=np.sqrt(2 * kB *temperature*gamma/dt)*np.dot(G_inv_t , eta)
       compute_msd=np.mean((r[i-1])**2)
       rumore[i]=noise
       r[i] = r[i-1] +rumore[i-1]*dt
       msd[i]=compute_msd

plt.figure(figsize=(10, 6))
# Grafico delle componenti della posizione rispetto al tempo
plt.subplot(2, 1, 1)
plt.plot(time, r[:,0], label='componente x')
plt.plot(time, r[:,1], label="componente y")
plt.xlabel('t')
plt.legend()
plt.ylabel('posizione')

plt.subplot(2, 1, 2)
plt.plot(time, msd)
plt.xlabel('t')
plt.ylabel('msd')
plt.show()



