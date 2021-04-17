# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 15:44:02 2018

@author: dell
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 13:47:56 2018

@author: dell
"""
from IPython import get_ipython
get_ipython().magic('reset -sf') 

import matplotlib.pyplot as plt
import numpy as np
import random
#import mpi4py 
import sys

from mpi4py import MPI
comm = MPI.COMM_WORLD
N_proc = comm.Get_size()
idx_proc = comm.Get_rank()
print (MPI.COMM_WORLD.size)
print('number of proc', N_proc)



L = 8
N_MC = 1000
comm = MPI.COMM_WORLD
N_proc = comm.Get_size()
idx_proc = comm.Get_rank()
N_warm = 3000
plt.figure(figsize=(8,8))
s = np.random.choice([1, -1], size=(L, L))
s_old = s
#plt.imshow(s)
#np.zeros((L,L))
J = 1.0
kB = 1.0
T = 8.0
NT = 8
ET = np.zeros(NT)
MT = np.zeros(NT)
CT = np.zeros(NT)
chiT = np.zeros(NT)

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


# begin main
Ea = np.zeros(N_MC)
Ma = np.zeros(N_MC)
Ca = np.zeros(N_MC)
chia = np.zeros(N_MC)
TT = np.zeros(NT)
T = 0.
for iT in range(0,NT):
    T = T+0.5
    TT[iT] = T
    print('T = ', T)
    acc = 0.
# warm up session
    for iter in range(N_warm):
        # begin sweep
        for sw in range(0,L*L):
            x = sw%L
            y = sw//L
            rdn = random.uniform(0, 1)
            prob = min([np.exp(-dE(s,L,x,y)/(kB*T)),1])
            if(  rdn < prob):#flip or not
                s[x,y] = -s[x,y]
                

    for iter in range(N_MC):
        # begin sweep
        for sw in range(0,L*L):
            x = sw%L
            y = sw//L
            rdn = random.uniform(0, 1)
            dffE = dE(s,L,x,y)
            if(abs(dffE -0)<1.0e-9):
                prob = random.uniform(0, 1)
            else:
                prob = min([np.exp(dffE/(kB*T)),1])
                #xxx = np.exp(dffE/(2*kB*T))
                #prob = xxx/(xxx + 1./xxx)
                
            if(  rdn < prob):#flip or not
                s[x,y] = -s[x,y]
                acc = acc+1./(N_MC*L*L)
        Ea[iter] = energy(s,L)
        Ma[iter] = abs(magnet(s,L))
        Ca[iter] = pow(Ea[iter],2)
        chia[iter] = Ma[iter]*Ma[iter]

        
#plt.imshow(s)
    Eave = sum(Ea)/N_MC
    Mave = sum(Ma)/N_MC
    Cave = pow(kB*T,-2)*( sum(Ca)/N_MC - pow(Eave,2))*L*L
    chiave = pow(kB*T,-1)*( sum(chia)/N_MC - pow(Mave,2) )*L*L
    print('average E = ', Eave)
    print('average M = ', Mave)
    print('average C = ', Cave)
    print('average chi = ', chiave)
    print('acceptance ratio = ', acc)
    ET[iT] = Eave
    MT[iT] = Mave
    CT[iT] = Cave
    chiT[iT] = chiave
    
plt.figure(figsize=(8,8))
plt.plot(TT, ET)
plt.xlabel('T', fontsize=20) 
plt.ylabel('Energy', fontsize=20) 
plt.figure(figsize=(8,8))
plt.plot(TT, MT)
plt.xlabel('T', fontsize=20) 
plt.ylabel('abs Mag', fontsize=20) 
plt.figure(figsize=(8,8))
plt.plot(TT, CT)
plt.xlabel('T', fontsize=20) 
plt.ylabel('Capacity', fontsize=20) 
plt.figure(figsize=(8,8))
plt.plot(TT, chiT)
plt.xlabel('T', fontsize=20) 
plt.ylabel('mag susceptibility', fontsize=20) 

print("size L = ",L, "max C = ",np.max(CT))
print("size L = ",L, "max at M = ",MT[4])


    
    

        
            




        