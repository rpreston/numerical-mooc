# -*- coding: utf-8 -*-
"""
Created on Sat Sep 27 19:44:27 2014

@author: Ryan
"""

import numpy
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams["font.family"] = "serif"
rcParams["font.size"] = 16

nx = 41
dx = 2./(nx-1)
nt = 10
dt = 0.001
c = 1.

u = numpy.ones(nx)
u[.5/dx : 1/dx + 1] = 2

#x = numpy.linspace(0,2,nx)

for n in range (1,nt):
    un = u.copy()
    u[1:] = un[1:]-un[1:]*dt/dx*(un[1:]-un[0:-1])
    u[0] = 1.0

plt.plot(numpy.linspace(0,2,nx), u, color='#003366', ls='--', lw=3)
plt.ylim(0,2.5);
'''
for n in range(1,nt):
    un = u.copy()
    for i in range (1, nx):
        u[i] = un[i]-c*dt/dx*(un[i]-un[i-1])

plt.plot(numpy.linspace(0,2,nx), u, color='#003366', ls='--', lw=3)
plt.ylim(0,2.5);
'''