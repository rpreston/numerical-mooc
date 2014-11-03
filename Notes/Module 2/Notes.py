# -*- coding: utf-8 -*-
"""
Created on Sun Nov 02 13:16:58 2014

@author: Ryan
"""

import numpy                       
import matplotlib.pyplot as plt    

from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def linearconv(nx):
    """Solve the linear convection equation.
    
    Solves the equation d_t u + c d_x u = 0 where 
    * the wavespeed c is set to 1
    * the domain is x \in [0, 2]
    * 20 timesteps are taken, with \Delta t computed using the CFL 0.5
    * the initial data is the hat function
    
    Produces a plot of the results
    
    Parameters
    ----------
    
    nx : integer
        number of internal grid points
        
    Returns
    -------
    
    None : none
    """
    dx = 2./(nx-1)
    nt = 20    
    c = 1
    sigma = .5
    
    dt = sigma*dx

    u = numpy.ones(nx) 
    u[.5/dx : 1/dx+1]=2

    un = numpy.ones(nx)

    for n in range(nt):  
        un = u.copy() 
        u[1:] = un[1:] -c*dt/dx*(un[1:] -un[0:-1]) 
        u[0] = 1.0
        
    plt.plot(numpy.linspace(0,2,nx), u, color='#003366', ls='--', lw=3)
    plt.ylim(0,2.5);