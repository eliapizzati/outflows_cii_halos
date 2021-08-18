# -*- coding: utf-8 -*-
"""
Created on Sun May  9 18:46:45 2021

@author: anna
"""



import matplotlib.pyplot as plt

import numpy as np

A_C = 2.69e-4 # carbon percentage for Z=1
A_H = 0.7 # hydrogen percentage for Z=1

Zeta = 1.

x_CII = 1.

x_e = 1.

T = np.logspace(1, 5, 1000)

n = 1.

epsilon1 = 3e-27 /1.4e-4 * (0.42*x_CII/1e-3) * np.exp(-92./T)

epsilon2 = 7.9e-20 * np.exp(-92./T)*T**(-0.5) # /  92**0.5 

epsilon3 = 0.67*T**0.13 *T**(-0.5)*8.6e-6/2 * np.exp(-92./T) * 0.0079*1.6e-12
fig, ax = plt.subplots()

ax.plot(T,  A_C / A_H * epsilon1)


ax.plot(T, A_C/A_H * epsilon2)

ax.plot(T, A_C/A_H * epsilon3)


ax.set_xscale("log")

ax.set_yscale("log")