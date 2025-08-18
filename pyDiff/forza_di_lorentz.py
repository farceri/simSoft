import numpy as np
import matplotlib.pyplot as plt
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

    
def verlet_integration_noise (mass, gamma, dt, N, numero_particelle, temperature, kB=1.0):
    r=np.zeros(shape=(N,2,numero_particelle))
    v=np.zeros(shape=(N,2,numero_particelle))
    msd=np.zeros(shape=(N,2,numero_particelle))
    msd_media=np.zeros(shape=(N,2))
    msd_xy_media=np.zeros(shape=(N,1))
    msd_yx_media=np.zeros(shape=(N,1))
    msd_xy = np.zeros(shape=(N, numero_particelle))
    msd_yx = np.zeros(shape=(N, numero_particelle))
    for i in range(1, N):
        #integrazione eulero
        for j in range(0, numero_particelle):
            noise=np.sqrt(2 * kB *temperature *gamma /dt)*np.array([np.random.normal(0,1),np.random.normal(0,1)])
            forza_trascinamento=-gamma * v[i-1,:,j]
            langevin_force=forza_trascinamento+noise
            v[i,:,j] = v[i-1,:,j] + 0.5*(langevin_force/ mass) * dt
            r[i,:,j] = r[i-1,:,j] + v[i-1,:,j]*dt +noise*dt
            msd[i,:,j]=(r[i,:,j]-r[0,:,j])**2
            msd_xy[i,j] = (r[i,0,j] - r[0,0,j])*(r[i,1,j] - r[0,1,j])
            msd_yx[i,j] = (r[i,1,j] - r[0,1,j])*(r[i,0,j] - r[0,0,j])
        # Media sulle particelle
        msd_media[i,:]=np.mean(msd[i,:,:], axis=1) 
        msd_xy_media[i] = np.mean(msd_xy[i,:], axis=0)
        msd_yx_media[i] = np.mean(msd_yx[i,:], axis=0)
    #
    return r, v, msd, msd_xy, msd_yx, msd_media, msd_xy_media, msd_yx_media


r, v, msd, msd_xy, msd_yx, msd_media, msd_xy_media, msd_yx_media = verlet_integration_noise(
    mass, gamma, dt, N, numero_particelle, temperature, kB)

def line_with_intercept_zero(x, b, q):
    return x*b + q


#Traiettoria
plt.plot(r[:,0,0], r[:,1,0], label='particella 1')
plt.plot(r[:,0,1], r[:,1,1], label='particella 2')
plt.plot(r[:,0,2], r[:,1,2], label='particella 3')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Percorso di una particella')
plt.axis('equal')
plt.legend()
plt.grid()
plt.show()


#CALCOLO MSD
plt.figure(figsize=(6, 6))
plt.subplot(2, 1, 1)
plt.plot(time, msd_media[:,0])
plt.xlabel('t')
plt.ylabel('msd_x')
plt.title('msd_x')
plt.subplot(2, 1, 2)
plt.xscale('log')
plt.yscale('log')
plt.plot(time, msd_media[:,0])
plt.xlabel('t')
plt.ylabel('msd_x')
plt.show()

plt.figure(figsize=(6, 6))
plt.subplot(2, 1, 1)
plt.plot(time, msd_media[:,1])
plt.xlabel('t')
plt.ylabel('msd_y')
plt.title('msd_y')
plt.subplot(2, 1, 2)
plt.xscale('log')
plt.yscale('log')
plt.plot(time,msd_media[:,1])
plt.xlabel('t')
plt.ylabel('msd_y')
plt.show()

plt.figure(figsize=(6, 6))
plt.subplot(2, 1, 1)
plt.plot(time, msd_xy_media[:,:])
plt.xlabel('t')
plt.ylabel('msd_xy')
plt.title('msd_xy,m_yx')
plt.subplot(2, 1, 2)
plt.plot(time, msd_yx_media[:,:])
plt.xlabel('t')
plt.ylabel('msd_yx')
plt.show()


#SELEZIONE VALORI MSD PER STIMA D
plt.figure(figsize=(6, 6))
plt.subplot(2, 1, 1)
plt.plot(time, (msd_media[:,0]/time))
plt.xlabel('t')
plt.ylabel('msd_xx/t')
plt.title('msd_xx/t,m_yy/t')
plt.subplot(2, 1, 2)
plt.plot(time,(msd_media[:,1]/time))
plt.xlabel('t')
plt.ylabel('msd_yy/t')
plt.show()

plt.figure(figsize=(6, 6))
plt.subplot(2, 1, 1)
plt.plot(time, (msd_xy_media[:,0]/time))
plt.xlabel('t')
plt.ylabel('msd_xy/t')
plt.title('msd_xy/t,m_yx/t')
plt.subplot(2, 1, 2)
plt.plot(time,(msd_yx_media[:,0]/time))
plt.xlabel('t')
plt.ylabel('msd_yx/t')
plt.show()



#STIMA DEL COEFFICIENTE DI DIFFUSIONE DAL FIT
#STIMA COEFFICIENTE DI DIFFUSIONE SECONDO METODO
tempo_D=700
D_xx=msd_media[tempo_D:N-1,0]/(2*(time[tempo_D:N-1]))
D_xx_media=np.mean(D_xx)
D_yy=msd_media[tempo_D:N-1,1]/(2*(time[tempo_D:N-1]))
D_yy_media=np.mean(D_yy)
D_xy=(msd_xy_media[tempo_D:N-1,0])/(2*(time[tempo_D:N-1]))
D_xy_media=np.mean(D_xy)
D_yx=(msd_yx_media[tempo_D:N-1,0])/(2*(time[tempo_D:N-1]))
D_yx_media=np.mean(D_yx)

#fit lineare componente xx con intercetta
parameters, error = curve_fit(line_with_intercept_zero, time[tempo_D:N-1], msd_media[tempo_D:N-1,0])
bx = parameters[0]
qx = parameters[1]
error_bx = np.sqrt(error[0,0])
error_qx = np.sqrt(error[1,1])

print(f"y = {bx:.2f} * x")
plt.yscale('log')
plt.xscale('log')
plt.plot(time[tempo_D:N-1], msd_media[tempo_D:N-1,0], color='blue', label='msd__xx' )
plt.plot(time[tempo_D:N-1], bx * time[tempo_D:N-1] + qx, color='purple', label='Fit y=a_x*t +q')
plt.plot(time[tempo_D:N-1], D_xx_media*2 * time[tempo_D:N-1], color='green', label='Fit y=D_x*2*t')
plt.legend()
plt.xlabel("t")
plt.ylabel("msd_xx")
plt.title("Fit lineare con intercetta (y = ax +q)")
plt.legend()
D_xx_fit=bx/2
error_D_xx_fit = error_bx/2
print('D_xx=',D_xx_media)
print('D_xx_fit_ +/- error =',D_xx_fit, ' +/- ', error_D_xx_fit)
print('q =',bx, ' +/- ', error_bx)
plt.show()

#fit lineare componente yy con intercetta
parameters, error = curve_fit(line_with_intercept_zero, time[tempo_D:N-1], msd_media[tempo_D:N-1,0])
bx = parameters[0]
qx = parameters[1]
error_bx = np.sqrt(error[0,0])
error_qx = np.sqrt(error[1,1])

print(f"y = {bx:.2f} * x")
plt.yscale('log')
plt.xscale('log')
plt.plot(time[tempo_D:N-1], msd_media[tempo_D:N-1,1], color='blue', label='msd__yy' )
plt.plot(time[tempo_D:N-1], bx * time[tempo_D:N-1] + qx, color='purple', label='Fit y=a_x*t +q')
plt.plot(time[tempo_D:N-1], D_yy_media*2 * time[tempo_D:N-1], color='green', label='Fit y=D_y*2*t')
plt.legend()
plt.xlabel("t")
plt.ylabel("msd_yy")
plt.title("Fit lineare con intercetta (y = ax +q)")
plt.legend()
#ricavo il coefficiente di diffusione dal fit
D_yy_fit=bx/2
error_D_yy_fit = error_bx/2
print('D_yy=',D_yy_media)
print('D_yy_fit_ +/- error =',D_yy_fit, ' +/- ', error_D_yy_fit)
print('q =',bx, ' +/- ', error_bx)
plt.show()

#fit lineare componente xy con intercetta
parameters, error = curve_fit(line_with_intercept_zero, time[tempo_D:N-1],msd_xy_media[tempo_D:N-1,0])
bx = parameters[0]
qx = parameters[1]
error_bx = np.sqrt(error[0,0])
error_qx = np.sqrt(error[1,1])

print(f"y = {bx:.2f} * x")
plt.yscale('log')
plt.xscale('log')
plt.plot(time[tempo_D:N-1], msd_xy_media[tempo_D:N-1], color='blue', label='msd__xy' )
plt.plot(time[tempo_D:N-1], bx * time[tempo_D:N-1] + qx, color='purple', label='Fit y=a_x*t +q')
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
print('q =',bx, ' +/- ', error_bx)
plt.show()

#fit lineare componente yx con intercetta
parameters, error = curve_fit(line_with_intercept_zero, time[tempo_D:N-1],msd_yx_media[tempo_D:N-1,0])
bx = parameters[0]
qx = parameters[1]
error_bx = np.sqrt(error[0,0])
error_qx = np.sqrt(error[1,1])

print(f"y = {bx:.2f} * x")
plt.yscale('log')
plt.xscale('log')
plt.plot(time[tempo_D:N-1], msd_yx_media[tempo_D:N-1], color='blue', label='msd__yx' )
plt.plot(time[tempo_D:N-1], bx * time[tempo_D:N-1] + qx, color='purple', label='Fit y=a_x*t +q')
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
print('q =',bx, ' +/- ', error_bx)
plt.show()
plt.close()

#VERIFICA COEFFICIENTE D DIFFUSIONE
D_verifica=(kB*temperature)/(gamma*mass)
D=np.sqrt((D_xx)**2+ (D_yy)**2)
print('D_verifica=',D_verifica)
print('D=',D)

#SALVATAGGIO COEFFICIENTI DI DIFFUSIONE
fname_no_lorentz = "C:\\Users\\viviana\\Desktop\\simulazioni\\tesi\\coefficenti_di_diffusione.out"
np.savetxt(fname_no_lorentz, (D_xx_media, D_yy_media, D_xy_media, D_yx_media, D_xx_fit,D_yy_fit,D_xy_fit,D_yx_fit,))



#FORZA DI LORENTZ
def verlet_integration_noise_fl (mass, gamma, dt, N, numero_particelle, temperature, Bz1, q, kB=1.0):
    r=np.zeros(shape=(N,2,numero_particelle))
    v=np.zeros(shape=(N,2,numero_particelle)) 
    msd=np.zeros(shape=(N,2,numero_particelle))
    msd_xy=np.zeros(shape=(N,numero_particelle))
    msd_yx=np.zeros(shape=(N,numero_particelle))
    msd_media=np.zeros(shape=(N,2))
    msd_xy_media=np.zeros(shape=(N,1))
    msd_yx_media=np.zeros(shape=(N,1))
    print("Verlet integration check, Bz = ", Bz1)
    for i in range(1, N):
        #integrazione eulero
        for j in range(0, numero_particelle):
            noise=np.sqrt(2 * kB * temperature * gamma / dt) * np.array(
                [np.random.normal(0,1), np.random.normal(0,1)]
            )
            G=np.array([[gamma,-q*Bz1],[q*Bz1,gamma]])
            forza_trascinamento=np.dot(-G , v[i-1,:,j])
            langevin_force=forza_trascinamento+noise
            v[i,:,j] = v[i-1,:,j] + 0.5*(langevin_force/ mass) * dt
            r[i,:,j] = r[i-1,:,j] + v[i-1,:,j]*dt +noise*dt
            msd[i,:,j]=(r[i,:,j] - r[0,:,j])**2
            msd_xy[i,j] = (r[i,0,j] - r[0,0,j])*(r[i,1,j] - r[0,1,j])
            msd_yx[i,j] = (r[i,1,j] - r[0,1,j])*(r[i,0,j] - r[0,0,j])
        # Media sulle particelle
        msd_media[i,:]=np.mean(msd[i,:,:], axis=1) 
        msd_xy_media[i] = np.mean(msd_xy[i], axis=0)
        msd_yx_media[i] = np.mean(msd_yx[i], axis=0)
        
    #
    return r, v, msd, msd_xy, msd_yx, msd_media, msd_xy_media, msd_yx_media

def line_with_intercept_zero_fl(x, b, q):
    return x*b + q

gamma=0.05 
t_diff=1/gamma
print('t_diff=',t_diff)
temperature=1e-3
kB=1.0
mass=1
sigma=np.sqrt(2*mass*gamma*kB*temperature)
numero_particelle=100
t_finale=800

q=1
Bz_values = np.array([0.00001,0.0001,0.001,0.01,0.1,0.5,0.8,1.0,1.2])
#print('Bz_values', Bz_values)
with open("C:\\Users\\viviana\\Desktop\\simulazioni\\tesi\\msd_data.npy", "wb") as f:
    np.save(f, Bz_values)
   
    
D_xx_media_list = []
D_yy_media_list = []
D_xy_media_list = []
D_yx_media_list = []
D_xx_fit_list = []
D_yy_fit_list = []
D_xy_fit_list = []
D_yx_fit_list = []

for k in range(len(Bz_values)):
    Bz = Bz_values[k]
    B = np.array([0,0,Bz])
    print('Bz=',Bz)
    N = 10000
    dt = 0.1  
    tempo_totale=N*dt           
    time = np.arange(0, N*dt, dt)

    r, v, msd, msd_xy, msd_yx, msd_media, msd_xy_media, msd_yx_media = verlet_integration_noise_fl(
         mass, gamma, dt, N, numero_particelle, temperature, Bz, q)
    
    # Save the data
    with open("C:\\Users\\viviana\\Desktop\\simulazioni\\tesi\\msd_data.npy", "ab") as f:
        np.save(f, time[:])
        np.save(f, msd_media[:,0])
        np.save(f, msd_media[:,1])
        np.save(f, msd_xy_media[:])
        np.save(f, msd_yx_media[:])
    
    
    
    #Traiettoria
    plt.plot(r[:,0,0], r[:,1,0])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Percorso di una particella, B_z={0:.5f}'.format(Bz))
    plt.axis('equal')
    plt.grid()
    plt.show()
    plt.close()

    #CALCOLO MSD
    plt.figure(figsize=(6, 6))
    plt.subplot(2, 1, 1)
    plt.plot(time, msd_media[:,0])
    plt.xlabel('t')
    plt.ylabel('msd_xx')
    plt.title('msd_xx,B_z={0:.5f}'.format(Bz))
    plt.subplot(2, 1, 2)
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(time,msd_media[:,0])
    plt.xlabel('t')
    plt.ylabel('msd_xx')
    plt.show()

    plt.figure(figsize=(6, 6))
    plt.subplot(2, 1, 1)
    plt.plot(time, msd_media[:,1])
    plt.xlabel('t')
    plt.ylabel('msd_yy')
    plt.title('msd_yy, B_z={0:.5f}'.format(Bz))
    plt.subplot(2, 1, 2)
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(time,msd_media[:,1])
    plt.xlabel('t')
    plt.ylabel('msd_yy')
    plt.show()

    plt.figure(figsize=(6, 6))
    plt.subplot(2, 1, 1)
    plt.plot(time, msd_xy_media[:,:])
    plt.xlabel('t')
    plt.ylabel('msd_xy')
    plt.title('msd_xy, msd_yx. B_z={0:.5f}'.format(Bz))
    plt.subplot(2, 1, 2)
    plt.plot(time, msd_yx_media[:,:])
    plt.xlabel('t')
    plt.ylabel('msd_yx')
    plt.show()

    #SELEZIONE VALORI MSD PER STIMA D
    plt.figure(figsize=(6, 6))
    plt.subplot(2, 1, 1)
    plt.plot(time, (msd_media[:,0]/time))
    plt.xlabel('t')
    plt.ylabel('msd_xx/t')
    plt.title('msd_xx/t, msd_yy/t, B_z={0:.5f}'.format(Bz))
    plt.subplot(2, 1, 2)
    plt.plot(time,(msd_media[:,1]/time))
    plt.xlabel('t')
    plt.ylabel('msd_yy/t')
    plt.show()

    plt.figure(figsize=(6, 6))
    plt.subplot(2, 1, 1)
    plt.plot(time, (msd_xy_media[:,0]/time))
    plt.xlabel('t')
    plt.ylabel('msd_xy/t')
    plt.title('msd_xy/t, msd_yx/t, B_z={0:.5f}'.format(Bz))
    plt.subplot(2, 1, 2)
    plt.plot(time,(msd_yx_media[:,0]/time))
    plt.xlabel('t')
    plt.ylabel('msd_yx/t')
    plt.show()
    

    #STIMA DEL COEFFICIENTE DI DIFFUSIONE DAL FIT
    #STIMA COEFFICIENTE DI DIFFUSIONE SECONDO METODO
    tempo_D=200
    tempo_f=800
    D_xx=msd_media[tempo_D:tempo_f,0]/(2*(time[tempo_D:tempo_f]))
    D_xx_media=np.mean(D_xx)
    D_xx_media_list.append(D_xx_media)
    D_yy=msd_media[tempo_D:tempo_f,1]/(2*(time[tempo_D:tempo_f]))
    D_yy_media=np.mean(D_yy)
    D_yy_media_list.append(D_yy_media)
    D_xy=(msd_xy_media[tempo_D:tempo_f,0]*msd_yx_media[tempo_D:tempo_f,0])/(2*(time[tempo_D:tempo_f]))
    D_xy_media=np.mean(D_xy)
    D_xy_media_list.append(D_xy_media)
    D_yx=(msd_yx_media[tempo_D:tempo_f,0]*msd_xy_media[tempo_D:tempo_f,0])/(2*(time[tempo_D:tempo_f]))
    D_yx_media=np.mean(D_yx)
    D_yx_media_list.append(D_yx_media)

    #fit lineare componente xx con intercetta
    
    parameters, error = curve_fit(line_with_intercept_zero_fl, time[tempo_D:tempo_f], msd_media[tempo_D:tempo_f,0])
    bx = parameters[0]
    qx = parameters[1]
    error_bx = np.sqrt(error[0,0])
    error_qx = np.sqrt(error[1,1])

    print(f"y = {bx:.5f} * x")
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(time[tempo_D:tempo_f], bx * time[tempo_D:tempo_f] + qx, color='purple', label='Fit y=a_x*t +q')
    plt.plot(time[tempo_D:tempo_f], D_xx_media*2 * time[tempo_D:tempo_f], color='green', label='Fit y=D_xx*2*t')
    plt.legend()
    plt.xlabel("t")
    plt.ylabel("msd_xx")
    plt.title("Fit lineare con intercetta (y = ax +q), B_z={0:.5f}".format(Bz))
    plt.legend()
    D_xx_fit=bx/2
    error_D_xx_fit = error_bx/2
    print('D_xx=',D_xx_media)
    print('D_xx_fit_ +/- error =',D_xx_fit, ' +/- ', error_D_xx_fit)
    D_xx_fit_list.append(D_xx_fit)
    print('q =',bx, ' +/- ', error_bx)
    plt.show()
    plt.close()

    parameters, error = curve_fit(line_with_intercept_zero_fl, time[tempo_D:tempo_f], msd_media[tempo_D:tempo_f,0])
    bx = parameters[0]
    qx = parameters[1]
    error_bx = np.sqrt(error[0,0])
    error_qx = np.sqrt(error[1,1])

    #fit lineare componente yy con intercetta
    print(f"y = {bx:.5f} * x")
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(time[tempo_D:tempo_f], bx * time[tempo_D:tempo_f] + qx, color='purple', label='Fit y=a_x*t +q')
    plt.plot(time[tempo_D:tempo_f], D_yy_media*2 * time[tempo_D:tempo_f], color='green', label='Fit y=D_yy*2*t')
    plt.legend()
    plt.xlabel("t")
    plt.ylabel("msd_yy")
    plt.title("Fit lineare con intercetta (y = ax +q),B_z={0:.5f}".format(Bz))
    plt.legend()
    #ricavo il coefficiente di diffusione dal fit
    D_yy_fit=bx/2
    error_D_yy_fit = error_bx/2
    print('D_yy=',D_yy_media)
    print('D_yy_fit_ +/- error =',D_yy_fit, ' +/- ', error_D_yy_fit)
    D_yy_fit_list.append(D_yy_fit)
    print('q =',bx, ' +/- ', error_bx)
    plt.show()

    #fit lineare componente xy con intercetta
    parameters, error = curve_fit(line_with_intercept_zero_fl, time[tempo_D:tempo_f],msd_xy_media[tempo_D:tempo_f,0])
    bx = parameters[0]
    qx = parameters[1]
    error_bx = np.sqrt(error[0,0])
    error_qx = np.sqrt(error[1,1])

    print(f"y = {bx:.5f} * x")
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(time[tempo_D:tempo_f], msd_xy_media[tempo_D:tempo_f],color='blue', label='dati' )
    plt.plot(time[tempo_D:tempo_f], bx * time[tempo_D:tempo_f] + qx, color='purple', label='Fit y=a_x*t +q')
    plt.plot(time[tempo_D:tempo_f], D_xy_media*2 * time[tempo_D:tempo_f], color='green', label='Fit y=D_xy*2*t')
    plt.legend()
    plt.xlabel("t")
    plt.ylabel("msd_xy")
    plt.title("Fit lineare con intercetta (y = ax +q),B_z={0:.5f}".format(Bz))
    plt.legend()
    #ricavo il coefficiente di diffusione dal fit
    D_xy_fit=bx/2
    error_D_xy_fit = error_bx/2
    print('D_xy=',D_xy_media)
    print('D_xy_fit_ +/- error =',D_xy_fit, ' +/- ', error_D_xy_fit)
    D_xy_fit_list.append(D_xy_fit)
    print('q =',bx, ' +/- ', error_bx)
    plt.show()


    #fit lineare componente yx con intercetta
    parameters, error = curve_fit(line_with_intercept_zero_fl, time[tempo_D:tempo_f],msd_yx_media[tempo_D:tempo_f,0])
    bx = parameters[0]
    qx = parameters[1]
    error_bx = np.sqrt(error[0,0])
    error_qx = np.sqrt(error[1,1])

    print(f"y = {bx:.5} * x")
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(time[tempo_D:tempo_f], bx * time[tempo_D:tempo_f] + qx, color='purple', label='Fit y=a_x*t +q')
    plt.plot(time[tempo_D:tempo_f], D_yx_media*2 * time[tempo_D:tempo_f], color='green', label='Fit y=D_yx*2*t')
    plt.legend()
    plt.xlabel("t")
    plt.ylabel("msd_yx")
    plt.title("Fit lineare con intercetta (y = ax +q), B_z={0:.5f}".format(Bz))
    plt.legend()
    D_yx_fit=bx/2
    error_D_yx_fit = error_bx/2
    print('D_yx=',D_yx_media)
    print('D_yx_fit_ +/- error =',D_yx_fit, ' +/- ', error_D_yx_fit)
    D_yx_fit_list.append(D_yx_fit)
    print('q =',bx, ' +/- ', error_bx)
    plt.show()
    plt.close()


#SALVATAGGIO COEFFICIENTI DI DIFFUSIONE
fname_lorentz = "C:\\Users\\viviana\\Desktop\\simulazioni\\tesi\\coefficenti_di_diffusione_fl.out"
np.savetxt(fname_lorentz, (Bz_values, D_xx_media_list, D_yy_media_list, D_xy_media_list, D_yx_media_list,
     D_xx_fit_list, D_yy_fit_list, D_xy_fit_list, D_yx_fit_list))


#COEFFICIENTE DI DIFFUSIONE RISPETTO AL CAMPO MAGNETICO
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.xscale('log')
plt.plot(Bz_values, D_xx_media_list, label='D')
plt.plot(Bz_values, D_xx_fit_list, label='D_fit')
plt.xlabel('B_z')
plt.ylabel('D_xx')
plt.legend()
plt.subplot(2, 1, 2)
plt.xscale('log')
plt.plot(Bz_values, D_yy_media_list, label='D')
plt.plot(Bz_values, D_yy_fit_list, label='D_fit')
plt.xlabel('B_z')
plt.ylabel('D_yy')
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.xscale('log')
#plt.plot(Bz_values, D_xy_media_list, label='D')
plt.plot(Bz_values, D_xy_fit_list, label='D_fit')
plt.xlabel('B_z')
plt.ylabel('D_xy')
plt.legend()
plt.subplot(2, 1, 2)
plt.xscale('log')
#plt.plot(Bz_values, D_yx_media_list, label='D')
plt.plot(Bz_values, D_yx_fit_list, label='D_fit')
plt.xlabel('B_z')
plt.ylabel('D_yx')
plt.legend()
plt.show()



