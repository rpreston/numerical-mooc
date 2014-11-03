# -*- coding: utf-8 -*-
"""
Created on Sun Nov 02 14:15:18 2014

@author: Ryan
"""

import numpy                       
import matplotlib.pyplot as plt    

from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

nx = 41
dx = 2./(nx-1)
nt = 20   
nu = 0.3   #the value of viscosity
sigma = .2 
dt = sigma*dx**2/nu 

x = numpy.linspace(0,2,nx)

u = numpy.ones(nx)      
u[.5/dx : 1/dx+1]=2  

un = numpy.ones(nx)

for n in range(nt):  
    un = u.copy() 
    u[1:-1] = un[1:-1] + nu*dt/dx**2*(un[2:] -2*un[1:-1] +un[0:-2]) 
        
plt.plot(numpy.linspace(0,2,nx), u, color='#003366', ls='--', lw=3)
plt.ylim(0,2.5);