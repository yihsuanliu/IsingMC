
import numpy as np
import random
#import mpi4py 
from mpi4py import MPI
comm = MPI.COMM_WORLD
N_proc = comm.Get_size()
idx_proc = comm.Get_rank()
print (MPI.COMM_WORLD.size)
#if(idx_proc==0):
print('number of proc', N_proc)


#define the parameters
L = 2
N_MC = 64000
N_MC_proc = N_MC//N_proc
N_warm = 8000
s = np.random.choice([1, -1], size=(L, L))
s_old=s
J = 1.0
kB = 1.0
N_T = 21 # number of temperature ticks
T = np.linspace(0.5,5,N_T)# temperature
Esave = np.zeros(N_T)
Msave = np.zeros(N_T)
Csave = np.zeros(N_T)
chisave= np.zeros(N_T)
#Hamiltonian

def energy(s,L):
    E = 0
    Energy = 0.0
    global J
    for i in range(L):
        ip1 = (i+1)%L
        for j in range(L):   
            jp1 = (j+1)%L
            E += s[i,j]*( s[ip1,j] +s[i,jp1])
    Energy = (-1)*E*J
    return Energy/(L*L)

def magnet(s, L):
    M = 0.0
    for i in range(L):        
        for j in range(L):   
            M += s[i,j]    
    return M/(L*L)
def dE(s,L,i,j):
    x1 = (i + 1)%L
    x0 = (i - 1 + L)%L
    y1 = (j + 1)%L
    y0 = (j - 1 + L)%L
    NB = s[x1,j] + s[x0,j] + s[i,y1] + s[i,y0]
    s0 = s[i,j]
    diff = -2*s0*NB*J
    #print(diff)
    return diff


# begin main program
for iT in range(0,N_T):
    iT2 = N_T-1-iT
    Temp = T[iT2]
    print('Temp = ', Temp)
    s = s_old
    acc = 0.
    E = 0.
    E2 = 0.
    M = 0
    M2 = 0.
# warm up session
    for iter in range(N_warm):
        # begin sweep
        for sw in range(0,L*L):
            x = sw%L
            y = sw//L
            rdn = random.uniform(0, 1)
            prob = min([np.exp(-dE(s,L,x,y)/(kB*Temp)),1])
            if(  rdn < prob):#flip or not
                s[x,y] = -s[x,y]
                

    for iter in range(N_MC_proc):
        # begin sweep
        for sw in range(0,L*L):
            x = sw%L
            y = sw//L
            rdn = random.uniform(0, 1)
            dffE = dE(s,L,x,y)
            if(abs(dffE -0)<1.0e-9):
                prob = random.uniform(0, 1)
            else:
                prob = min([np.exp(dffE/(kB*Temp)),1])
                #xxx = np.exp(dffE/(2*kB*T))
                #prob = xxx/(xxx + 1./xxx)
                
            if(  rdn < prob):#flip or not
                s[x,y] = -s[x,y]
                acc = acc+1./(N_MC*L*L)
                
        e = energy(s,L)
        E += e
        E2 += e*e
        m = abs(magnet(s,L))
        M += m
        M2 += m*m
        s_old = s
        
    #collect from each workder
    dummy = np.array([E], dtype ='d')
    dummy2 = np.empty(1, dtype ='d')
    comm.Allreduce(dummy,dummy2)
    E = dummy2
    dummy = np.array([E2], dtype ='d')
    dummy2 = np.empty(1, dtype ='d')
    comm.Allreduce(dummy,dummy2)
    E2 = dummy2
    dummy = np.array([M], dtype ='d')
    dummy2 = np.empty(1, dtype ='d')
    comm.Allreduce(dummy,dummy2)
    M = dummy2
    dummy = np.array([M2], dtype ='d')
    dummy2 = np.empty(1, dtype ='d')
    comm.Allreduce(dummy,dummy2)
    M2 = dummy2

    # calc the for temperatuer 
    Eave = E/N_MC
    Mave = M/N_MC
    Cave = pow(kB*Temp,-2)*( (E2)/N_MC - pow(Eave,2))*L*L
    chiave = pow(kB*Temp,-1)*( (M2)/N_MC - pow(Mave,2) )*L*L
    if(idx_proc==0):
        print('average E = ', Eave)
        print('average M = ', Mave)
        print('average C = ', Cave)
        print('average chi = ', chiave)
        print('acceptance ratio = ', acc)
    Esave[iT2] = Eave
    Msave[iT2] = Mave
    Csave[iT2] = Cave
    chisave[iT2] = chiave
np.savetxt('L'+str(L)+'E.dat',Esave)
np.savetxt('L'+str(L)+'M.dat',Msave)
np.savetxt('L'+str(L)+'C.dat',Csave)
np.savetxt('L'+str(L)+'chi.dat',chisave)
#%%
'''
import matplotlib.pyplot as plt
plt.figure(figsize=(8,8))
plt.plot(T, Esave)
plt.xlabel('T', fontsize=20) 
plt.ylabel('Energy', fontsize=20) 
plt.figure(figsize=(8,8))
plt.plot(T, Msave)
plt.xlabel('T', fontsize=20) 
plt.ylabel('abs Mag', fontsize=20) 
plt.figure(figsize=(8,8))
plt.plot(T, Csave)
plt.xlabel('T', fontsize=20) 
plt.ylabel('Capacity', fontsize=20) 
plt.figure(figsize=(8,8))
plt.plot(T, chisave)
plt.xlabel('T', fontsize=20) 
plt.ylabel('mag susceptibility', fontsize=20) 

#print("size L = ",L, "max C = ",np.max(CT))
#print("size L = ",L, "max at M = ",MT[4])
''' 
