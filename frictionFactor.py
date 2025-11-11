# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 08:42:05 2025

@author: frekih
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from flowProperties import q, d, mu, rho

f0 = 1
f1 = 101
df = 0.1

u = q / (3600 * np.pi * (d/2)**2)

N_sim = int((f1 - f0) / df) + 1

t_array = np.zeros(N_sim)
f_array = np.zeros(N_sim)

for i in np.arange(0.001, 0.05, 0.001):
    epsilon = i

r = epsilon / d # roughness [-]
v = mu / rho # kinematic viscosity [m2/s]
Re = u * d / v

for k in range(0, N_sim):
    
    dt = df * k
    
    def f(x):
        if Re < 4000:
            x = 64 / Re * dt
        else:
            return (-2 * np.log10((r / 3.7) + (2.51 / (Re * np.sqrt(x)))) - (1 / np.sqrt(x)))
        
    sol = fsolve(f, 0.001)
    ff = sol * dt
        
    f_array[k] = ff
    t_array[k] = dt
    
plt.Figure(figsize=[16, 12])
plt.plot(t_array, f_array, 'y', label='ff')
plt.xlim(0, 1)
plt.ylim(0.4, 0)
plt.grid()
plt.show()