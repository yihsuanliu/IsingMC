import matplotlib.pyplot as plt
import numpy as np
N_T = 21
T = np.linspace(0.5,5,N_T)# temperature
Lm = np.array([2, 4, 8, 16, 20])

#%% plot Energy
plt.figure(figsize=(8,8))
for iL in range(5):
    L = Lm[iL]
    Esave = np.genfromtxt('L'+str(L) + 'E.dat')
    plt.plot(T, Esave, label = 'L = ' +str(L))
    plt.xlabel('T', fontsize=20) 
    plt.ylabel('Energy', fontsize=20) 
plt.legend(loc='upper left', fontsize=20)     
#%% plot Magnetization
plt.figure(figsize=(8,8))
for iL in range(5):
    L = Lm[iL]
    Msave = np.genfromtxt('L'+str(L) + 'M.dat')
    plt.plot(T, Msave, label = 'L = ' +str(L))
    plt.xlabel('T', fontsize=20) 
    plt.ylabel('abs Mag', fontsize=20) 
plt.legend(loc='lower left', fontsize=20)     

# %% plot heat capacity
plt.figure(figsize=(8,8))
for iL in range(5):
    L = Lm[iL]
    Csave = np.genfromtxt('L'+str(L) + 'C.dat')
    plt.plot(T, Csave, label = 'L = ' +str(L))
    plt.xlabel('T', fontsize=20)
    plt.ylabel('Capacity', fontsize=20) 
plt.legend(loc='upper right', fontsize=20)     
plt.ylim(0,5)


#%% plot mag susceptibility
plt.figure(figsize=(8,8))
for iL in range(4):
    L = Lm[iL]
    chisave = np.genfromtxt('L'+str(L) + 'chi.dat')    
    plt.plot(T, chisave, label = 'L = ' +str(L))
    plt.ylabel('mag susceptibility', fontsize=20) 
plt.legend(loc='upper right', fontsize=20)     
plt.ylim(0,10)