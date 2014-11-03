# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 21:08:45 2014

@author: Ryan
"""
import numpy
import matplotlib.pyplot as plt

from matplotlib import animation

Vmax = 80. #km/hr
#Vmax = 136 #2nd part of questions
L = 11. #km
rhomax = 250. #cars/km
nx = 51
dt = 0.001 #hours or 0.06 min

dx = L/nx
x = numpy.linspace(0,L,nx)
rho0 = numpy.ones(nx)*10.
#rho0 = numpy.ones(nx)*20. #2nd part of questions
rho0[10:20] = 50.
nt = 1


def carVelocity(rho):
    return (Vmax*(1-rho/rhomax))

def trafficFlux(rho,Velocity):
    return (Velocity*rho)


def evaluateRho(t):
    #line.set_data(x,rho0)
    for n in range(1,int(t/dt)):
        rho = rho0.copy()
        rho0[1:] = rho[1:]-(dt/dx)*Vmax*(1.-(rho[1:])*(2./rhomax))*(rho[1:]-rho[0:-1])
        
'''  
def plotGraph(t):
    nt = int(t/dt)
    print nt
    fig = plt.figure(figsize=(8,5))
    ax = plt.axes(xlim=(0,55),ylim=(0,55))
    line = ax.plot([],[],color='#003366', ls='--', lw=3)[0]
   
    animation.FuncAnimation(fig,evaluateRho,frames= nt,interval=100)

'''