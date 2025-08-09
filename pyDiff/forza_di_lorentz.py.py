import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.stats import norm
import sys
import sympy as sp
from scipy.signal import butter, lfilter
from scipy.optimize import curve_fit


gamma=0.05 
t_diff=1/gamma
print('t_diff=',t_diff)
temperature=1e-3
kB=1.0
mass=1
sigma=np.sqrt(2*mass*gamma*kB*temperature)
numero_particelle=100


#INTEGRALE DELL'ENERGIA con il rumore(bidimensionale)
N = 10000  
dt = 0.1  
tempo_totale=N*dt           
time = np.arange(0, N*dt, dt)
v0=np.sqrt(kB * temperature /1)
v_0=np.array([v0,v0])
r_0=np.array([0,0])
r=np.zeros(shape=(N,2,numero_particelle))
v=np.zeros(shape=(N,2,numero_particelle))
E = np.empty(shape=(N,1,numero_particelle))  
msd=np.empty(shape=(N,2,numero_particelle))
msd_media=np.empty(shape=(N,2))
msd_xy_media=np.empty(shape=(N,1))
msd_yx_media=np.empty(shape=(N,1))

        
def verlet_integration_noise (r, v, E,msd,msd_media,msd_xy_media,msd_yx_media,mass, gamma,dt,N,temperature):
        n_part = np.shape(r)[2]
        for j in range(1, n_part):
            
            #integrazione verlet
            for i in range(1, N):
                noise=np.sqrt(2 * kB *temperature *gamma /dt)*np.array([np.random.normal(0,1),np.random.normal(0,1)])
                forza_trascinamento=-gamma * v[i-1,:,j]
                langevin_force=forza_trascinamento+noise
                #kinetic_energy = 0.5 *mass *((np.linalg.norm(v[i-1,:,j]))** 2)
                v[i,:,j] = v[i-1,:,j] + 0.5*(langevin_force/ mass) * dt
                r[i,:,j] = r[i-1,:,j] + v[i-1,:,j]*dt +noise*dt
                #E[i,:,j] =kinetic_energy
                msd[i,:,:]=((r[i,:,:]-r[0,:,:]))**2
                msd_media[i,:]=np.mean(msd[i,:,:], axis=1)
                msd_xy = (r[i,0,:] - r[0,0,:]) * (r[i,1,:] - r[0,1,:])
                msd_xy_media[i] = np.mean(msd_xy)
                msd_yx = (r[i,0,:] - r[0,0,:]) * (r[i,1,:] - r[0,1,:])
                msd_yx_media[i] = np.mean(msd_yx)

            
        return r, v, E,msd, msd_media,msd_xy_media, msd_yx_media
r[:,:,:], v[:,:,:], E[:,:,:],msd, msd_media,msd_xy_media,msd_yx_media,= verlet_integration_noise (r, v, E,msd, msd_media,msd_xy_media,msd_yx_media, mass, gamma,dt,N,temperature)


#Traiettoria
plt.plot(r[:,0,1], r[:,1,1],label='particella 1')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Percorso di una particella (Moto Browniano 2D)')
plt.legend()
plt.axis('equal')
plt.grid()
plt.show()


#CALCOLO MSD
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(time, msd_media[:,0])
plt.xlabel('t')
plt.ylabel('msd_xx')
plt.subplot(2, 1, 2)
plt.plot(np.log10(time), np.log10(msd_media[:,0]))
plt.xlabel('log(t)')
plt.ylabel('log(msd_xx)')
plt.show()

plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(time, msd_media[:,1])
plt.xlabel('t')
plt.ylabel('msd_yy')
plt.subplot(2, 1, 2)
plt.plot(np.log10(time), np.log10(msd_media[:,1]))
plt.xlabel('log(t)')
plt.ylabel('log(msd_yy)')
plt.show()

plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(time, msd_xy_media[:,:])
plt.xlabel('t')
plt.ylabel('msd_xy')
plt.subplot(2, 1, 2)
plt.plot(np.log10(time), np.log10(msd_xy_media[:,:]))
plt.xlabel('log(t)')
plt.ylabel('log(msd_xy)')
plt.show()

plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(time, msd_yx_media[:,:])
plt.xlabel('t')
plt.ylabel('msd_yx')
plt.subplot(2, 1, 2)
plt.plot(np.log10(time), np.log10(msd_yx_media[:,:]))
plt.xlabel('log(t)')
plt.ylabel('log(msd_yx)')
plt.show()


#SELEZIONE VALORI MSD PER STIMA D
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(time, (msd_media[:,0]/time))
plt.xlabel('t')
plt.ylabel('msd_xx/t')
plt.subplot(2, 1, 2)
plt.plot(time,(msd_media[:,1]/time))
plt.xlabel('t')
plt.ylabel('msd_yy/t')
plt.show()

plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(time, (msd_xy_media[:,0]/time))
plt.xlabel('t')
plt.ylabel('msd_xy/t')
plt.subplot(2, 1, 2)
plt.plot(time,(msd_yx_media[:,0]/time))
plt.xlabel('t')
plt.ylabel('msd_yx/t')
plt.show()



#STIMA DEL COEFFICIENTE DI DIFFUSIONE DAL FIT
#STIMA COEFFICIENTE DI DIFFUSIONE SECONDO METODO
tempo_D=400
D_xx=msd_media[tempo_D:N-1,0]/(2*tempo_totale)
D_xx_media=np.mean(D_xx)
D_yy=msd_media[tempo_D:N-1,1]/(2*tempo_totale)
D_yy_media=np.mean(D_yy)
D_xy=(msd_xy_media[tempo_D:N-1,0])/(2*tempo_totale)
D_xy_media=np.mean(D_xy)
D_yx=(msd_yx_media[tempo_D:N-1,0])/(2*tempo_totale)
D_yx_media=np.mean(D_yx)


#fit lineare componente xx con intercetta
def line_with_intercept_zero(x, b, q):
      return x*b + q

parameters, error = curve_fit(line_with_intercept_zero, time[tempo_D:N-1], msd_media[tempo_D:N-1,0])
bx = parameters[0]
q = parameters[1]
error_bx = np.sqrt(error[0,0])
error_q = np.sqrt(error[1,1])

print(f"y = {bx:.2f} * x")
plt.yscale('log')
plt.xscale('log')
plt.plot(time[tempo_D:N-1], bx * time[tempo_D:N-1] + q, color='purple', label='Fit y=a_x*t +q')
plt.plot(time[tempo_D:N-1], D_xx_media*2 * time[tempo_D:N-1], color='green', label='Fit y=D_xx*2*t')
plt.legend()
plt.xlabel("t")
plt.ylabel("msd_x")
plt.title("Fit lineare con intercetta (y = ax +q)")
plt.legend()
D_xx_fit=bx/2
error_D_xx_fit = error_bx/2
print('D_xx=',D_xx_media)
print('D_xx_fit_ +/- error =',D_xx_fit, ' +/- ', error_D_xx_fit)
plt.show()



#fit lineare componente yy con intercetta
def line_with_intercept_zero(x, b, q):
      return x*b + q

parameters, error = curve_fit(line_with_intercept_zero, time[tempo_D:N-1], msd_media[tempo_D:N-1,0])
bx = parameters[0]
q = parameters[1]
error_bx = np.sqrt(error[0,0])
error_q = np.sqrt(error[1,1])

print(f"y = {bx:.2f} * x")
plt.yscale('log')
plt.xscale('log')
plt.plot(time[tempo_D:N-1], bx * time[tempo_D:N-1] + q, color='purple', label='Fit y=a_x*t +q')
plt.plot(time[tempo_D:N-1], D_yy_media*2 * time[tempo_D:N-1], color='green', label='Fit y=D_yy*2*t')
plt.legend()
plt.xlabel("t")
plt.ylabel("msd_y")
plt.title("Fit lineare con intercetta (y = ax +q)")
plt.legend()
#ricavo il coefficiente di diffusione dal fit
D_yy_fit=bx/2
error_D_yy_fit = error_bx/2
print('D_yy=',D_yy_media)
print('D_yy_fit_ +/- error =',D_yy_fit, ' +/- ', error_D_yy_fit)
plt.show()


#fit lineare componente xy con intercetta
def line_with_intercept_zero(x, b, q):
      return x*b + q

parameters, error = curve_fit(line_with_intercept_zero, time[tempo_D:N-1],msd_xy_media[tempo_D:N-1,0])
bx = parameters[0]
q = parameters[1]
error_bx = np.sqrt(error[0,0])
error_q = np.sqrt(error[1,1])

print(f"y = {bx:.2f} * x")
plt.yscale('log')
plt.xscale('log')
plt.plot(time[tempo_D:N-1], bx * time[tempo_D:N-1] + q, color='purple', label='Fit y=a_x*t +q')
plt.plot(time[tempo_D:N-1], D_xy_media*2 * time[tempo_D:N-1], color='green', label='Fit y=D_xy*2*t')
plt.legend()
plt.xlabel("t")
plt.ylabel("msd_xy")
plt.title("Fit lineare con intercetta (y = ax +q)")
plt.legend()
#ricavo il coefficiente di diffusione dal fit
D_xy_fit=bx/2
error_D_xy_fit = error_bx/2
print('D_xy=',D_xy_media)
print('D_xy_fit_ +/- error =',D_xy_fit, ' +/- ', error_D_xy_fit)
plt.show()



#fit lineare componente yx con intercetta
def line_with_intercept_zero(x, b, q):
      return x*b + q

parameters, error = curve_fit(line_with_intercept_zero, time[tempo_D:N-1],msd_yx_media[tempo_D:N-1,0])
bx = parameters[0]
q = parameters[1]
error_bx = np.sqrt(error[0,0])
error_q = np.sqrt(error[1,1])

print(f"y = {bx:.2f} * x")
plt.yscale('log')
plt.xscale('log')
plt.plot(time[tempo_D:N-1], bx * time[tempo_D:N-1] + q, color='purple', label='Fit y=a_x*t +q')
plt.plot(time[tempo_D:N-1], D_yx_media*2 * time[tempo_D:N-1], color='green', label='Fit y=D_yx*2*t')
plt.legend()
plt.xlabel("t")
plt.ylabel("msd_yx")
plt.title("Fit lineare con intercetta (y = ax +q)")
plt.legend()
D_yx_fit=bx/2
error_D_yx_fit = error_bx/2
print('D_yx=',D_yx_media)
print('D_yx_fit_ +/- error =',D_yx_fit, ' +/- ', error_D_yx_fit)
plt.show()

#VERIFICA COEFFICIENTE D DIFFUSIONE
D_verifica=(kB*temperature)/(gamma*mass)
print('D_verifica=',D_verifica)





#FORZA DI LORENTZ

def verlet_integration_noise (r, v, E,msd,msd_media,msd_xy_media,msd_yx_media,mass, gamma,dt,N,temperature):
      n_part = np.shape(r)[2]
      for j in range(1, n_part):
            #integrazione verlet
            for i in range(1, N):
                  noise=np.sqrt(2 * kB *temperature *gamma /dt)*np.array([np.random.normal(0,1),np.random.normal(0,1)])
                  G=np.array([[gamma,-q*B[2]],[q*B[2],gamma]])
                  forza_trascinamento=np.dot(-G , v[i-1,:,j])
                  langevin_force=forza_trascinamento+noise
                  kinetic_energy = 0.5 *mass *((np.linalg.norm(v[i-1,:,j]))** 2)
                  v[i,:,j] = v[i-1,:,j] + 0.5*(langevin_force/ mass) * dt
                  r[i,:,j] = r[i-1,:,j] + v[i-1,:,j]*dt +noise*dt
                  E[i,:,j] =kinetic_energy
                  msd[i,:,:]=(r[i,:,:]-r[0,:,:])**2
                  msd_media[i,:]=np.mean(msd[i,:,:], axis=1) 
                  msd_xy = (r[i,0,:] - r[0,0,:]) * (r[i,1,:] - r[0,1,:])
                  msd_xy_media[i] = np.mean(msd_xy)
                  msd_yx = (r[i,0,:] - r[0,0,:]) * (r[i,1,:] - r[0,1,:])
                  msd_yx_media[i] = np.mean(msd_yx)
      
      return r, v, E,msd,msd_media,msd_xy_media,msd_yx_media
      


q=1
v_fl=np.array([1,2,0])
B_z0=0.4
B_0=np.array([0,0,B_z0])
B_max=1.2
Bz_values = np.arange(B_z0, B_max, 0.2)
print('Bz_values', Bz_values)
D_xx_media_list = []
D_yy_media_list = []
D_xy_media_list = []
D_yx_media_list = []
D_xx_fit_list = []
D_yy_fit_list = []
D_xy_fit_list = []
D_yx_fit_list = []
for Bz in Bz_values:
      B = np.array([0,0,Bz])
      print('Bz=',Bz)
      forza_di_lorentz=np.cross(q*v_fl,B)
      print('La forza di LORENTZ è:',forza_di_lorentz)




      #INTEGRALE DELL'ENERGIA con forza di lorentz(bidimensionale)
      N = 10000   
      dt = 0.1  
      tempo_totale=N*dt           
      time = np.arange(0, N*dt, dt)
      v0=np.sqrt(kB *temperature /1)
      v_0=np.array([v0,v0])
      r_0=np.array([0,0])
      r=np.zeros(shape=(N,2,numero_particelle))
      v=np.zeros(shape=(N,2,numero_particelle))
      E = np.empty(shape=(N,1,numero_particelle))  
      r_2=np.empty(shape=(N,2,numero_particelle))
      msd=np.empty(shape=(N,2,numero_particelle))
      msd_media=np.empty(shape=(N,2))
      msd_xy_media=np.empty(shape=(N,1))
      msd_yx_media=np.empty(shape=(N,1))

      r[:,:,:], v[:,:,:], E[:,:,:],msd,msd_media,msd_xy_media,msd_yx_media= verlet_integration_noise (r, v, E,msd, msd_media,msd_xy_media,msd_yx_media, mass, gamma,dt,N,temperature)

            
      


      #Traiettoria
      plt.plot(r[:,0,1], r[:,1,1],label='particella 1')
      plt.xlabel('x')
      plt.ylabel('y')
      plt.title('Percorso di una particella (Moto Browniano 2D con forza di Lorentz), B_z={0:.2f}'.format(Bz))
      plt.legend()
      plt.axis('equal')
      plt.grid()
      plt.show()


      #CALCOLO MSD
      plt.figure(figsize=(10, 6))
      plt.subplot(2, 1, 1)
      plt.plot(time, msd_media[:,0])
      plt.xlabel('t')
      plt.ylabel('msd_xx')
      plt.title('msd_xx, B_z={0:.2f}'.format(Bz))
      plt.subplot(2, 1, 2)
      plt.plot(np.log10(time), np.log10(msd_media[:,0]))
      plt.xlabel('log(t)')
      plt.ylabel('log(msd_xx)')
      plt.show()

      plt.figure(figsize=(10, 6))
      plt.subplot(2, 1, 1)
      plt.plot(time, msd_media[:,1])
      plt.xlabel('t')
      plt.ylabel('msd_yy')
      plt.title('msd_yy, B_z={0:.2f}'.format(Bz))
      plt.subplot(2, 1, 2)
      plt.plot(np.log10(time), np.log10(msd_media[:,1]))
      plt.xlabel('log(t)')
      plt.ylabel('log(msd_yy)')
      plt.show()

      plt.figure(figsize=(10, 6))
      plt.subplot(2, 1, 1)
      plt.plot(time, msd_xy_media[:,:])
      plt.xlabel('t')
      plt.ylabel('msd_xy')
      plt.title('msd_xy, B_z={0:.2f}'.format(Bz))
      plt.subplot(2, 1, 2)
      plt.plot(np.log10(time), np.log10(msd_xy_media[:,:]))
      plt.xlabel('log(t)')
      plt.ylabel('log(msd_xy)')
      plt.show()

      plt.figure(figsize=(10, 6))
      plt.subplot(2, 1, 1)
      plt.plot(time, msd_yx_media[:,:])
      plt.xlabel('t')
      plt.ylabel('msd_yx')
      plt.title('msd_yx, B_z={0:.2f}'.format(Bz))
      plt.subplot(2, 1, 2)
      plt.plot(np.log10(time), np.log10(msd_yx_media[:,:]))
      plt.xlabel('log(t)')
      plt.ylabel('log(msd_yx)')
      plt.show()


      #SELEZIONE VALORI MSD PER STIMA D
      plt.figure(figsize=(10, 6))
      plt.subplot(2, 1, 1)
      plt.plot(time, (msd_media[:,0]/time))
      plt.xlabel('t')
      plt.ylabel('msd_xx/t')
      plt.title('msd_xx/t, msd_yy/t, B_z={0:.2f}'.format(Bz))
      plt.subplot(2, 1, 2)
      plt.plot(time,(msd_media[:,1]/time))
      plt.xlabel('t')
      plt.ylabel('msd_yy/t')
      plt.show()

      plt.figure(figsize=(10, 6))
      plt.subplot(2, 1, 1)
      plt.plot(time, (msd_xy_media[:,0]/time))
      plt.xlabel('t')
      plt.ylabel('msd_xy/t')
      plt.title('msd_xy/t, msd_yx/t, B_z={0:.2f}'.format(Bz))
      plt.subplot(2, 1, 2)
      plt.plot(time,(msd_yx_media[:,0]/time))
      plt.xlabel('t')
      plt.ylabel('msd_yx/t')
      plt.show()
      



      #STIMA DEL COEFFICIENTE DI DIFFUSIONE DAL FIT
      #STIMA COEFFICIENTE DI DIFFUSIONE SECONDO METODO
      D_xx=msd_media[tempo_D:N-1,0]/(2*tempo_totale)
      D_xx_media=np.mean(D_xx)
      D_xx_media_list.append(D_xx_media)
      D_yy=msd_media[tempo_D:N-1,1]/(2*tempo_totale)
      D_yy_media=np.mean(D_yy)
      D_yy_media_list.append(D_yy_media)
      D_xy=(msd_xy_media[tempo_D:N-1,0]*msd_yx_media[tempo_D:N-1,0])/(2*tempo_totale)
      D_xy_media=np.mean(D_xy)
      D_xy_media_list.append(D_xy_media)
      D_yx=(msd_yx_media[tempo_D:N-1,0]*msd_xy_media[tempo_D:N-1,0])/(2*tempo_totale)
      D_yx_media=np.mean(D_yx)
      D_yx_media_list.append(D_yx_media)

      #fit lineare componente xx con intercetta
      def line_with_intercept_zero(x, b, q):
            return x*b + q

      parameters, error = curve_fit(line_with_intercept_zero, time[tempo_D:N-1], msd_media[tempo_D:N-1,0])
      bx = parameters[0]
      q = parameters[1]
      error_bx = np.sqrt(error[0,0])
      error_q = np.sqrt(error[1,1])

      print(f"y = {bx:.2f} * x")
      plt.yscale('log')
      plt.xscale('log')
      plt.plot(time[tempo_D:N-1], bx * time[tempo_D:N-1] + q, color='purple', label='Fit y=a_x*t +q')
      plt.plot(time[tempo_D:N-1], D_xx_media*2 * time[tempo_D:N-1], color='green', label='Fit y=D_xx*2*t')
      plt.legend()
      plt.xlabel("t")
      plt.ylabel("msd_x")
      plt.title("Fit lineare con intercetta (y = ax +q), B_z={0:.2f}".format(Bz))
      plt.legend()
      D_xx_fit=bx/2
      error_D_xx_fit = error_bx/2
      print('D_xx=',D_xx_media)
      print('D_xx_fit_ +/- error =',D_xx_fit, ' +/- ', error_D_xx_fit)
      D_xx_fit_list.append(D_xx_fit)
      plt.show()



      #fit lineare componente yy con intercetta
      def line_with_intercept_zero(x, b, q):
            return x*b + q

      parameters, error = curve_fit(line_with_intercept_zero, time[tempo_D:N-1], msd_media[tempo_D:N-1,0])
      bx = parameters[0]
      q = parameters[1]
      error_bx = np.sqrt(error[0,0])
      error_q = np.sqrt(error[1,1])

      print(f"y = {bx:.2f} * x")
      plt.yscale('log')
      plt.xscale('log')
      plt.plot(time[tempo_D:N-1], bx * time[tempo_D:N-1] + q, color='purple', label='Fit y=a_x*t +q')
      plt.plot(time[tempo_D:N-1], D_yy_media*2 * time[tempo_D:N-1], color='green', label='Fit y=D_yy*2*t')
      plt.legend()
      plt.xlabel("t")
      plt.ylabel("msd_y")
      plt.title("Fit lineare con intercetta (y = ax +q),B_z={0:.2f}".format(Bz))
      plt.legend()
      #ricavo il coefficiente di diffusione dal fit
      D_yy_fit=bx/2
      error_D_yy_fit = error_bx/2
      print('D_yy=',D_yy_media)
      print('D_yy_fit_ +/- error =',D_yy_fit, ' +/- ', error_D_yy_fit)
      D_yy_fit_list.append(D_yy_fit)
      plt.show()


      #fit lineare componente xy con intercetta
      def line_with_intercept_zero(x, b, q):
            return x*b + q

      parameters, error = curve_fit(line_with_intercept_zero, time[tempo_D:N-1],msd_xy_media[tempo_D:N-1,0])
      bx = parameters[0]
      q = parameters[1]
      error_bx = np.sqrt(error[0,0])
      error_q = np.sqrt(error[1,1])

      print(f"y = {bx:.2f} * x")
      plt.yscale('log')
      plt.xscale('log')
      plt.plot(time[tempo_D:N-1], bx * time[tempo_D:N-1] + q, color='purple', label='Fit y=a_x*t +q')
      plt.plot(time[tempo_D:N-1], D_xy_media*2 * time[tempo_D:N-1], color='green', label='Fit y=D_xy*2*t')
      plt.legend()
      plt.xlabel("t")
      plt.ylabel("msd_xy")
      plt.title("Fit lineare con intercetta (y = ax +q),B_z={0:.2f}".format(Bz))
      plt.legend()
      #ricavo il coefficiente di diffusione dal fit
      D_xy_fit=bx/2
      error_D_xy_fit = error_bx/2
      print('D_xy=',D_xy_media)
      print('D_xy_fit_ +/- error =',D_xy_fit, ' +/- ', error_D_xy_fit)
      D_xy_fit_list.append(D_yx_fit)
      plt.show()



      #fit lineare componente yx con intercetta
      def line_with_intercept_zero(x, b, q):
            return x*b + q

      parameters, error = curve_fit(line_with_intercept_zero, time[tempo_D:N-1],msd_yx_media[tempo_D:N-1,0])
      bx = parameters[0]
      q = parameters[1]
      error_bx = np.sqrt(error[0,0])
      error_q = np.sqrt(error[1,1])

      print(f"y = {bx:.2f} * x")
      plt.yscale('log')
      plt.xscale('log')
      plt.plot(time[tempo_D:N-1], bx * time[tempo_D:N-1] + q, color='purple', label='Fit y=a_x*t +q')
      plt.plot(time[tempo_D:N-1], D_yx_media*2 * time[tempo_D:N-1], color='green', label='Fit y=D_yx*2*t')
      plt.legend()
      plt.xlabel("t")
      plt.ylabel("msd_yx")
      plt.title("Fit lineare con intercetta (y = ax +q), B_z={0:.2f}".format(Bz))
      plt.legend()
      D_yx_fit=bx/2
      error_D_yx_fit = error_bx/2
      print('D_yx=',D_yx_media)
      print('D_yx_fit_ +/- error =',D_yx_fit, ' +/- ', error_D_yx_fit)
      D_yx_fit_list.append(D_yx_fit)
      plt.show()


#COEFFICIENTE DI DIFFUSIONE RISPETTO AL CAMPO MAGNeTICO
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(Bz_values, D_xx_media_list, label='D')
plt.plot(Bz_values, D_xx_fit_list, label='D_fit')
plt.xlabel('B_z')
plt.ylabel('D_xx')
plt.subplot(2, 1, 2)
plt.plot(Bz_values, D_yy_media_list, label='D')
plt.plot(Bz_values, D_yy_fit_list, label='D_fit')
plt.xlabel('B_z')
plt.ylabel('D_yy')
plt.show()

plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(Bz_values, D_xy_media_list, label='D')
plt.plot(Bz_values, D_xy_fit_list, label='D_fit')
plt.xlabel('B_z')
plt.ylabel('D_xy')
plt.subplot(2, 1, 2)
plt.plot(Bz_values, D_xy_media_list, label='D')
plt.plot(Bz_values, D_xy_fit_list, label='D_fit')
plt.xlabel('B_z')
plt.ylabel('D_xy')
plt.show()
