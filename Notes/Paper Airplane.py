# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 20:31:29 2014

@author: Ryan
"""
from math import sin, cos, log, ceil
import numpy
import matplotlib.pyplot as plt


from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

# model parameters:
g = 9.8      # gravity in m s^{-2}
v_t = 4.9   # trim velocity in m s^{-1}   
C_D = 1/5.  # drag coefficient --- or D/L if C_L=1
C_L = 1.0    # for convenience, use C_L = 1

### set initial conditions ###
v0 = v_t     # start at the trim velocity (or add a delta)
theta0 = 0.0 # initial angle of trajectory
x0 = 0.0     # horizotal position is arbitrary
y0 = 1000.0  # initial altitude



def f(u):
    """Returns the right-hand side of the phugoid system of equations.
    
    Parameters
    ----------
    u : array of float
        array containing the solution at time n.
        
    Returns
    -------
    dudt : array of float
        array containing the RHS given u.
    """
    
    v = u[0]
    theta = u[1]
    x = u[2]
    y = u[3]
    return numpy.array([-g*sin(theta) - C_D/C_L*g/v_t**2*v**2,
                      -g*cos(theta)/v + g/v_t**2*v,
                      v*cos(theta),
                      v*sin(theta)])
                      
def euler_step(u, f, dt):
    """Returns the solution at the next time-step using Euler's method.
    
    Parameters
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the right hand-side of the system of equation.
    dt : float
        time-increment.
    
    Returns
    -------
    u_n_plus_1 : array of float
        approximate solution at the next time step.
    """
    print u+dt*f(u)
    return u + dt * f(u)
    
T = 100.0                          # final time
dt = 0.1                           # time increment
N = int(T/dt) + 1                  # number of time-steps
t = numpy.linspace(0.0, T, N)      # time discretization

# initialize the array containing the solution for each time-step
u = numpy.empty((N, 4))
# fill 1st element with initial values


def nominalRun():
    
    theta0 = 0.0
    u[0] = numpy.array([v_t, theta0, x0, y0])
    
    
    for n in range(N-1):
        u[n+1] = euler_step(u[n],f,dt)
    x=u[:,2]
    y=u[:,3]
    plt.figure(figsize=(8,6))
    plt.grid(True)
    plt.xlabel(r'x', fontsize=18)
    plt.ylabel(r'y', fontsize=18)
    plt.title('Glider Nominal trajectory, flight time = %.2f' % T, fontsize=18)
    plt.plot(x,y,"k-",lw=2)
    plt.plot(t,x,'b-',lw=2)
    plt.plot(t,y,'r',lw=2)

def monteCarlo(iterations,maxArg,minArg):
    xm = numpy.empty((iterations,N))
    ym = numpy.empty((iterations,N))
    
    v0=v_t
    theta0 = minArg
    
    
    for i in range(0,iterations):
        theta0 += (maxArg-minArg)/(iterations)
        u = numpy.empty((N,4))
        u[0] = numpy.array([v0, theta0, x0, y0])
        
        for n in range(N-1):
            
            u[n+1] = euler_step(u[n], f, dt)
        
        x = u[:,2]
        y = u[:,3]
        for p in range(N-1):
            xm[i][p] = x[p]
            ym[i][p] = y[p]
    plt.figure(figsize=(8,6))
    plt.grid(True)
    plt.xlabel(r'x', fontsize=18)
    plt.ylabel(r'y', fontsize=18)
    plt.title('Glider trajectory, flight time = %.2f' % T, fontsize=18)
    for temp in range(0,iterations):
        plt.plot(xm[temp][:N-1],ym[temp][:N-1], 'k-',lw=1);
        plt.plot(t[:N-1],xm[temp][:N-1],'b',lw=1)
        plt.plot(t[:N-1],ym[temp][:N-1],'r',lw=1)
    
monteCarlo(50,numpy.pi,0.0)

        

# get the glider's position with respect to the time
#monteCarlo(50,7,2)
# visualization of the path
