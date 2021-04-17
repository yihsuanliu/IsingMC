# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 14:31:27 2018

@author: dell
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress



data = np.array([[10,1.116,2.576],
        [16,1.118,5.417],
        [20,1.36269,6.807],
        [24,1.284,7.72]])


x = np.log(data[:,0])

logchi = np.log(data[:,2])

plt.figure(figsize=(8,8))
plt.plot(x, logchi)
plt.xlabel('log(L)', fontsize=20) 
plt.ylabel('log(chi)', fontsize=20) 
    

