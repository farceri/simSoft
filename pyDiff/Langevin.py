import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from numba import njit



def langevin(tau, c, x0, dt, tstop):
    rng = np.random.default_rng()
    mu = np.exp(-(1 / tau) * dt)
    sigma = np.sqrt((c * tau / 2) * (1 - np.exp(-(2 / tau) * dt)))

    num_steps = int(tstop / dt)
    X = np.zeros(num_steps)
    T = np.zeros(num_steps)
    Y = np.zeros(num_steps)
    X[0] = x0
    T[0] = 0
    Y[0] = 0

    # Vectorized Langevin equation
    N = rng.normal(0, 1, size=num_steps - 1)  # Generate all random numbers at once
    for i in range(1, num_steps):
        X[i] = X[i - 1] * mu + sigma * N[i - 1]
        T[i] = T[i - 1] + dt

    # Vectorized calculations for med, var, devst, medy, vary
    i = np.arange(num_steps)
    med = x0 * np.exp(-(i * dt) / tau)
    var = (c * tau / 2) * (1 - np.exp(-2 * (i * dt) / tau))
    devst = np.sqrt(var)

    # Vectorized calculation for Y
    Y[1:] = np.cumsum(X[:-1] * dt)

    medy = x0 * tau * (1 - np.exp(-(i * dt) / tau))
    vary = c * tau**2 * (i * dt - 2 * tau * (1 - np.exp(-(i * dt) / tau)) + (tau / 2) * (1 - np.exp(-2 * (i * dt) / tau)))

    return T, X, med, devst, Y, medy, vary


'''
def langevin_direct(tau, c, x0, t_indices):
    """
    Calculates the Langevin equation solution directly at specified time steps.
    """
    rng = np.random.default_rng()
    T = np.array(t_indices)  # Time points at which to calculate the solution

    # Calculate X directly at the given time points
    X = rng.normal(x0*np.exp(-(T)/tau),(c*tau/2)*(1-np.exp(-2*(T)/tau)))

    return T, X




#@numba.jit("Tuple((f8[:],f8[:],f8[:],f8[:],f8[:],f8[:],f8[:]))(i8,i8,i8,i8,i8)",nopython=True,nogil=True)
def langevin(tau,c,x0,dt,tstop):
    rng = np.random.default_rng()
    mu = np.exp(-(1 / tau) * dt)
    sigma = np.sqrt((c * tau / 2) * (1 - np.exp(-(2 / tau) * dt)))

    X = np.zeros(int(tstop / dt))
    T = np.zeros(int(tstop / dt))
    Y=np.zeros(int(tstop/dt))
    X[0] = x0
    T[0] = 0
    Y[0]=0

    for i in range(1, int(tstop / dt)):
        N = rng.normal(0, 1)
        X[i] = X[i - 1] * mu + sigma * N
        T[i] = T[i - 1] + dt

    med = np.zeros(int(tstop / dt))
    for i in range(len(med)):
        med[i] = x0 * np.exp(-(i * dt) / tau)

    var = np.zeros(int(tstop / dt))
    for i in range(len(var)):
        var[i] = (c * tau / 2) * (1 - np.exp(-2 * (i * dt) / tau))

    devst = np.sqrt(var)

    for i in range(len(Y)):
        Y[i]=Y[i-1]+X[i-1]*dt

    medy = np.zeros(int(tstop / dt))
    for i in range(len(Y)):
        medy[i]=x0*tau*(1-np.exp(-(i*dt)/tau))

    vary=np.zeros(int(tstop / dt))
    for i in range(len(Y)):
        vary[i]=c*tau**2*(i*dt-2*tau*(1-np.exp(-(i*dt)/tau))+(tau/2)*(1-np.exp(-2*(i*dt)/tau)))


    return T,X,med,devst,Y,medy,vary
'''

'''
def potenza(tau,c):
    nu=np.logspace(-4,6,100)
    S=[]
    for i in range(len(nu)):
        S.append((2*c*tau**2)/(1+(2*np.pi*tau*nu[i])**2))
    return nu,S
'''



def potenza(tau, c):
    nu = np.logspace(-4, 6, 100)
    S = (2 * c * tau**2) / (1 + (2 * np.pi * tau * nu)**2)  # Vectorized calculation
    return nu, S




def main():



    tau=float(input("Dammi tau: \t"))
    c = float(input("Dammi c: \t"))
    x0 = float(input("Dammi x0: \t"))
    dt = float(input("Dammi dt: \t"))
    tstop=float(input("Dammi tstop: \t"))
    yn=int(input("Traiettorie [1], Gaussiane METODO LENTO[0], Gaussiane METODO VELOCE[2]: \t"))

    #yn=2

    if(yn==1):
        N = int(input("Quante copie vuoi: \t"))
        figure, axis=plt.subplots(2,2)
        plt.rcParams.update({'font.size': 5})
        for i in range(N):
            axis[0,0].plot(langevin(tau, c, x0, dt, tstop)[0], langevin(tau, c, x0, dt, tstop)[1])
            axis[1,0].plot(langevin(tau, c, x0, dt, tstop)[0], langevin(tau, c, x0, dt, tstop)[4])
            axis[0, 1].plot(potenza(tau, c)[0], potenza(tau, c)[1])


        axis[0,0].plot(langevin(tau, c, x0, dt, tstop)[0], langevin(tau, c, x0, dt, tstop)[2], color='k')
        axis[0,0].plot(langevin(tau, c, x0, dt, tstop)[0], langevin(tau, c, x0, dt, tstop)[3], ls='--', color='r')
        axis[0,0].plot(langevin(tau, c, x0, dt, tstop)[0], -langevin(tau, c, x0, dt, tstop)[3], ls='--', color='r')
        axis[0,0].set_title("Traiettorie O-U 1D")

        axis[1,0].plot(langevin(tau, c, x0, dt, tstop)[0], langevin(tau, c, x0, dt, tstop)[5], color='k')
        axis[1,0].plot(langevin(tau, c, x0, dt, tstop)[0],  langevin(tau, c, x0, dt, tstop)[5]+langevin(tau, c, x0, dt, tstop)[6], ls='--', color='r')
        axis[1,0].plot(langevin(tau, c, x0, dt, tstop)[0], langevin(tau, c, x0, dt, tstop)[5] -langevin(tau, c, x0, dt, tstop)[6], ls='--', color='r')
        axis[1,0].set_title("Integrale O-U 1D")

        axis[0,0].sharex(axis[1,0])



        axis[0, 1].set_title("Potenza dissipata all'eq")
        axis[0,1].set_yscale('log')
        axis[0,1].set_xscale('log')

        '''
        for i in range(2):
            for j in range(2):
                axis[i,j].grid(True, which="both")
        '''
        #axis[:,:].grid(True, which="both")



        plt.show()


    elif(yn==0):

       # R = int(input("Quante realizzazioni vuoi: \t"))
        R=int(tstop/dt)
        Xg3 = np.zeros(R)
        Xg23 = np.zeros(R)
        med3=0
        med23=0
        s3=0
        s23=0


        for i in range(R):
            Xg3[i] = langevin(tau, c, x0, dt, tstop)[1][int(tstop / 3)]
            Xg23[i] = langevin(tau, c, x0, dt, tstop)[1][int(2 * tstop / 3)]


        for i in range(R):
            med3+=Xg3[i]

            med23+=Xg23[i]
        med3=med3/R
        med23=med23/R
        for i in range(R):
            s3+=(Xg3[i]-med3)**2
            s23+=(Xg23[i]-med23)**2


        s3=s3/(R-1)
        s23=s23/(R-1)


        xaxes = ['Pos1', 'Pos2']
        yaxes = ['Prob', 'Prob']
        titles = ['Gaussiane a tstop/3 e 2tstop/3', '']
        data=[Xg3,Xg23]
        x_axis = np.arange(-np.max(data[:][0]), np.max(data[:][0]), 0.001)
        f, a = plt.subplots(2, 1)


        errm3=100*(x0*np.exp(-(tstop/3)/tau)-med3)/(x0*np.exp(-(tstop/3)/tau))
        errm23=100*(x0*np.exp(-(2*tstop/3)/tau)-med23)/(x0*np.exp(-(2*tstop/3)/tau))
        errv3=np.abs(100*((c*tau/2)*(1-np.exp(-2*(tstop/3)/tau))-s3)/((c*tau/2)*(1-np.exp(-2*(tstop/3)/tau))))
        errv23=np.abs(100*((c*tau/2)*(1-np.exp(-2*(2*tstop/3)/tau))-s23)/((c*tau/2)*(1-np.exp(-2*(2*tstop/3)/tau))))

        f=open("Par_gauss.txt","w")
        f.write("t~="+str(int(tstop/3))+"\n medth="+str(x0*np.exp(-(tstop/3)/tau))+"\n varth="+str(round((c*tau/2)*(1-np.exp(-2*(tstop/3)/tau)),10))+"\n med="+str(round(med3,4))+"\n var="+str(round(s3,4)))
        #f.write("\n errm="+str(100*(x0*np.exp(-(tstop/3)/tau))-med3))/(x0*np.exp(-(tstop/3)/tau))+"\n errv="+str(100*(((c*tau/2)*(1-np.exp(-2*(tstop/3)/tau)))-(s3))/(((c*tau/2)*(1-np.exp(-2*(tstop/3)/tau)))))
        f.write("\n errm="+str(errm3)+"\n errv="+str(errv3))
        f.write("\t \t \t")
        f.write("t~="+str(int(2*tstop/3))+"\n medth="+str(x0*np.exp(-(2*tstop/3)/tau))+"\n varth="+str(round((c*tau/2)*(1-np.exp(-2*(2*tstop/3)/tau)),10))+"\n med="+str(round(med23,4))+"\n var="+str(round(s23,4)))
        #f.write("\n errm="+str(100*(x0*np.exp(-(2*tstop/3)/tau))-(med23))/(x0*np.exp(-(2*tstop/3)/tau))+"\n errv="+str(100*(((c*tau/2)*(1-np.exp(-2*(2*tstop/3)/tau)))-(s23))/(((c*tau/2)*(1-np.exp(-2*(2*tstop/3)/tau))))))
        f.write("\n errm=" + str(errm23) + "\n errv=" + str(errv23))
        a = a.ravel()
        for idx, ax in enumerate(a):
            ax.hist(data[idx],weights=np.zeros_like(data[idx]) + 1. / data[idx].size)
            ax.set_title(titles[idx])
            ax.set_xlabel(xaxes[idx])
            ax.set_ylabel(yaxes[idx])


        bw1=(np.max(data[:][0])-np.min(data[:][0]))/10
        bw2 = (np.max(data[:][1]) - np.min(data[:][1]) )/ 10
        #moltiplico per larghezza bin*numero di conteggi totali dispense di doro
        a[0].plot(x_axis, bw1*norm.pdf(x_axis, x0*np.exp(-(tstop/3)/tau), (c*tau/2)*(1-np.exp(-2*(tstop/3)/tau))))
        a[1].plot(x_axis, bw2*norm.pdf(x_axis, x0 * np.exp(-(2*tstop / 3) / tau),(c * tau / 2) * (1 - np.exp(-2 * (2*tstop / 3) / tau))))
        #plt.tight_layout()
        plt.show()

    '''
    elif(yn==2):

        num_samples = int(input("Quanti sample vuoi?"))
        t_indices=input("A quali tempi (frazioni di tstop) vuoi valutare? [mettili separati da uno spazio]:").split()
        t_indices=np.array(t_indices,dtype=int)
        t_indices[:]=tstop/t_indices[:]

        #t_indices = [tstop / 3]  # Time points at which to calculate the solution

        # CALCOLO LE I PUNTI CHE MI INTERESSANO
        X_samples = np.zeros((num_samples, len(t_indices)))
        for i in range(num_samples):
            _, X = langevin_direct(tau, c, x0, t_indices)  # Call the direct calculation function
            X_samples[i, :] = X


        fig,axis=plt.subplots(1,len(t_indices),sharey=True)
        axis[0].set_ylabel("Prob")

        for i in range(len(t_indices)):
            bw= (np.max(X_samples[:,i]) - np.min(X_samples[:,i])) / 10
            x_axis = np.arange(-np.max(X_samples[:,i]), np.max(X_samples[:,i]), 0.001)
            axis[i].hist(X_samples[:,i],weights=np.zeros_like(X_samples[:,i]) + 1. / X_samples[:,i].size)
            axis[i].plot(x_axis,bw*norm.pdf(x_axis,x0*np.exp(-(t_indices[i])/tau), (c*tau/2)*(1-np.exp(-2*(t_indices[i])/tau))))
            axis[i].set_title("Tempo: "+str(t_indices[i]))
            axis[i].set_xlabel("Posizione")

        media=np.zeros(len(t_indices))
        var = np.zeros(len(t_indices))
        for i in range(len(t_indices)):
            media[i]= np.mean(X_samples[:,i])
            var[i] = np.var(X_samples[:,i])




        #print("media: " + str(media) + "\t devst: " + str(np.sqrt(var)))

        #OUTPUT SU FILE
        f=open("Parametri.txt","w")

        for i in range(len(media)):
            f.write("Media ["+str(i)+"]= "+str(media[i]))
            f.write("\t \t")
            f.write("Var ["+str(i)+"]= "+ str(var[i]))
            f.write("\t \t")
            f.write("Devst [" + str(i) + "]= " + str(np.sqrt(var[i])))
            f.write("\n \n")


        plt.show()
    '''



if __name__ == "__main__":
    main()

