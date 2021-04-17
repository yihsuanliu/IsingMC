# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 13:47:56 2018

@author: dell
"""

import matplotlib.pyplot as plt
import numpy as np
import random

L = 8
N_MC = 5000
N_warm = 2000
plt.figure(figsize=(8,8))
s = np.random.choice([1, -1], size=(L, L))
plt.imshow(s)
#np.zeros((L,L))
J = 1.0
kB = 1.0
T = 1
NT = 8
ET = np.array(NT)
MT = np.array(NT)
CT = np.array(NT)
chiT = np.array(NT)

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
        prob = min([np.exp(dE(s,L,x,y)/(kB*T)),1])
        if(  rdn < prob):#flip or not
            s[x,y] = -s[x,y]
    Ea[iter] = energy(s,L)
    Ma[iter] = abs(magnet(s,L))
    Ca[iter] = pow(Ea[iter],2)
    chia[iter] = Ma[iter]*Ma[iter]

        
#plt.imshow(s)
Eave = sum(Ea)/N_MC
Mave = sum(Ma)/N_MC
Cave = pow(kB*T,-2)*( sum(Ca)/N_MC - pow(Eave,2))
chiave = pow(kB*T,-1)*( sum(chia)/N_MC - pow(Mave,2) )
print('average E = ', Eave)
print('average M = ', Mave)
print('average C = ', Cave)
print('average chi = ', chiave)


plt.imshow(s)

        
            




        