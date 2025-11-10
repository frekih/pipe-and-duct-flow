# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 08:42:05 2025

@author: frekih
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from flowProperties import u, d, mu, rho

t0 = 0
t1 = 100
dt = 0.1

for i in np.arange(0.001, 0.05, 0.001):
    epsilon = i

r = epsilon / d # roughness [-]
v = mu / rho # kinematic viscosity [m2/s]
Re = u * d / v

def f(x):
    if Re < 4000:
        x = 64 / Re # for laminar flow
    else:
        return (-2 * np.log10((epsilon / 3.7 * d) + (2.51 / Re * np.sqrt(x))) - (1 / np.sqrt(x)))
    
sol = fsolve(f, 0.001)
print(sol)
    
plt.Figure(figsize=[16, 12])