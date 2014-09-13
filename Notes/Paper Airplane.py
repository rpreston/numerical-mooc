# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 20:31:29 2014

@author: Ryan
"""
from math import sin, cos, log, ceil, floor, sqrt
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
thetaInit = 0.0
theta0 = 0.0 # initial angle of trajectory
x0 = 0.0     # horizotal position is arbitrary
y0 = 10.0  # initial altitude


T = 11.                          # final time
dt = 0.1                           # time increment
N = int(T/dt) + 1                  # number of time-steps
t = numpy.linspace(0.0, T, N)      # time discretization

global_V0Max = 0.0
global_ThetaMax = 0.0

# initialize the array containing the solution for each time-step
u = numpy.empty((N, 4))

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
    if(u[3] <= 0.0):
        y = 0.0
        u[3] = 0.0
        u[2] = 0.0
        x = 0.0
    
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
    #print u+dt*f(u)
    if u[3]<= 0.0:
        return u
    return u + dt * f(u)
    

# fill 1st element with initial values

def paramRun(thetaArg,v0Arg):
    
    u[0] = numpy.array([v0Arg, thetaArg, x0, y0])
    
    n=0
    #while u[n][3] >= 0.4:
    for n in range(0,N-1):
        if u[n+1,3]<0.5 and u[n+1,3]>u[n,3]:
            u[n+1][3] = 0
            u[n+1][2] = u[n][2]
        else:
            u[n+1] = euler_step(u[n],f,dt)
        #n+=1
    
    for n in range(1,N):
        if u[n-1,3]<0.0:
            u[n,3]=0.0
            u[n,2]=u[n-1,2]
    
    x=u[:,2]
    y=u[:,3]
   
    plt.plot(x,y,lw=2)
    
def nominalRun():
    
    u[0] = numpy.array([v_t, thetaInit, x0, y0])
    
    n=0
    #while u[n][3] >= 0.4:
    for n in range(0,N-1):
        if u[n+1][3]<0.5 and u[n+1][3]>u[n][3]:
            u[n+1][3] = 0.0
            u[n+1][2] = u[n][2]
        else:
            u[n+1] = euler_step(u[n],f,dt) 
        #n+=1
    for n in range(1,N-1):
        if u[n-1,3]<0.0:
            u[n,3]=0.0
            u[n,2]=u[n-1,2]
        
    x=u[:,2]
    y=u[:,3]
   
    plt.plot(x,y,"b",lw=2)
    #plt.plot(t,x,'b-',lw=2)
    #plt.plot(t,y,'r',lw=2)

def monteCarlo(iterations,maxArgTheta,minArgTheta,maxArgV0,minArgV0):
    global xm
    global ym
    global global_V0Max
    global global_ThetaMax
    #global vm
    #global thetam
    floorSqrtIterations = floor(sqrt(iterations))
    
    xm = numpy.zeros((int(floorSqrtIterations**2),N))
    ym = numpy.zeros((int(floorSqrtIterations**2),N))
    vm = numpy.zeros((int(floorSqrtIterations**2),N))
    thetam = numpy.zeros((int(floorSqrtIterations**2),N))
    
    maxDistance = numpy.zeros((int(floorSqrtIterations**2)))
    maxTheta = numpy.zeros((int(floorSqrtIterations**2)))
    maxV0 = numpy.zeros((int(floorSqrtIterations**2)))
    maxArg = numpy.zeros((int(floorSqrtIterations**2)))
    #get floor sqrt iterations    
    
    global v0
    global theta0
    
    
    
    v0 = minArgV0
    theta0 = minArgTheta
    plt.figure(figsize=(8,8))
    for j in range(0,int(floorSqrtIterations)):
        v0+= (maxArgV0-minArgV0)/(floorSqrtIterations)
        theta0 = minArgTheta
        for i in range(0,int(floorSqrtIterations)):
            theta0 += (maxArgTheta-minArgTheta)/(floorSqrtIterations)
            u = numpy.zeros((N,4))
            u[0] = numpy.array([v0, theta0, x0, y0])
            
            paramRun(theta0,v0)
            nominalRun()
            #print 'u'
            #print u[0]
            #print 'i'
            #print i*j
            n = 0
            #while u[n][3] >= 0.0:
            for n in range(0,N-1):
                if u[n+1][3]<0.5 and u[n+1][3]>u[n][3]:
                    u[n+1][3] = 0.0
                    u[n+1][2] = u[n][2]
                else:
                    u[n+1] = euler_step(u[n],f,dt)
                #n+=1

            for n in range(1,N-1):
                if u[n-1,3]<0.0:
                    u[n,3]=0.0
                    u[n,2]=u[n-1,2]            
            
            x = u[:,2]
            y = u[:,3]
            v = u[:,0]
            the = u[:,1]
            
            for p in range(N-1):
                if y[p] <= 0.0:
                    x[p] = 0.0
                    y[p] = 0.0
            
            for p in range(N-1):
                xm[j*(int(floorSqrtIterations))+i][p] = x[p]
                ym[j*(int(floorSqrtIterations))+i][p] = y[p]
                vm[j*(int(floorSqrtIterations))+i][p] = v[p]
                thetam[j*(int(floorSqrtIterations))+i][p] = the[p]
            counter = 0
            for p in ym[j*(int(floorSqrtIterations))+i]:
                if p < 0.0:
                    xm[j*(int(floorSqrtIterations))+i][counter] = -1
                    ym[j*(int(floorSqrtIterations))+i][counter] = -1
                counter+=1
           
           
           
    print "# of Iterations"
    print floorSqrtIterations
    
    plt.figure(figsize=(8,8))
    plt.grid(True)
    plt.xlabel(r'x', fontsize=18)
    plt.ylabel(r'y', fontsize=18)
    plt.title('Glider trajectory, flight time = %.2f' % T, fontsize=18)
    #print floorSqrtIterations**2.
    for temp in range(0,int(floorSqrtIterations**2)):
        plt.plot(xm[temp][:N-1],ym[temp][:N-1], 'k-',lw=1);
    
    nominalRun()
    #plt.figure(figsize=(8,8))
    for temp in range(0,int(floorSqrtIterations**2)):
        maxDistance[temp] = max(xm[temp])
        maxArg[temp] = numpy.argmax(xm[temp])
        maxTheta[temp] = max(thetam[0])
        maxV0[temp] = max(vm[0])
        
    
        #plt.plot(t[:N-1],xm[temp][:N-1], 'k-',lw=1);
        
        #print temp
        #plt.plot(t[:N-1],xm[temp][:N-1],'b',lw=1)
        #plt.plot(t[:N-1],ym[temp][:N-1],'r',lw=1)
    print "Max X:"
    print max(maxDistance)
    print numpy.argmax(maxDistance)
    print maxDistance[numpy.argmax(maxDistance)]
    tempArg = numpy.argmax(maxDistance)
    
    #print "init V0 for max X"
    global_V0Max = maxV0[tempArg]
    #print "init theta for max X"
    global_ThetaMax = maxTheta[tempArg]
    
monteCarlo(5000.,thetaInit+0.8,thetaInit-0.8,v_t+0.8,v_t-0.8)


plt.figure(figsize=(8,8))
nominalRun()
paramRun(global_ThetaMax,global_V0Max)

# get the glider's position with respect to the time
#monteCarlo(50,7,2)
# visualization of the path
