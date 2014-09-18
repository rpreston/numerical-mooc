# -*- coding: utf-8 -*-
"""
Created on Tue Sep 09 20:18:55 2014

@author: Ryan
"""
import numpy
import matplotlib.pyplot as plt

#currentH = 0
h0 = 0
v0 = 0
totalRocketMass = 50 #kg
g = 9.81
rho = 1.091
r = 0.5
maxCrossSection = numpy.pi*r**2
t0 = 0.0
ve = 325
Cd = 0.15
mpo = 100
dt = 0.01
N = int(40.0/dt)

u = numpy.zeros((N,3))

#u = numpy.array([h0,v0])

def stepFunction(currentTime,startTime,endTime,highValue,lowValue):
    if currentTime>=startTime and currentTime<endTime :
        return highValue
    else:
        return lowValue


def getPropellantMass(currentTime):
    if currentTime<=5 :
        returnValue = (mpo-currentTime*stepFunction(currentTime,0.,5.,20.,0.))
    else:
        returnValue = 0
    return returnValue

def rocketEqn(argArray):
    mdotp = stepFunction(argArray[2],0.,5.,20.,0.)
    mp = getPropellantMass(argArray[2])
    numerator = mdotp*ve - (1./2.)*rho*argArray[1]*abs(argArray[1])*maxCrossSection*Cd
    denominator = totalRocketMass + mp
    return numpy.array([argArray[1],(-g) + numerator/denominator,1.])

def euler_Step(argArray, f, dt):
    return argArray + dt * f(argArray)

def getPosVel():
    global u
    count = 1
    u[0] = numpy.array([h0,v0,t0])
    while u[count-1][0]>0.0 or u[count-1][2]< 1.0:
        u[count] = euler_Step(u[count-1],rocketEqn,dt)
        count += 1
    for i in range(count,N):
        u[i][0] = u[i-1][0]
        u[i][1] = u[i-1][1]
        u[i][2] = u[i-1][2]
        
def plotAll():
    plt.figure(figsize=(8,8))
    plt.grid(True)
    plt.xlabel(r"Time", fontsize=18)
    plt.ylabel(r"Height", fontsize=18)
    plt.title("Rocket Trajectory", fontsize=18)
    plt.plot(u[:,2],u[:,0],u[:,2],u[:,1])

def getMaxAltData():
    print "Max Altitude Data"
    print "Altitude:", max(u[:,0])
    print "Velocity:",u[numpy.argmax(u[:,0]),1]
    print "Time:",u[numpy.argmax(u[:,0]),2]
    

def getMaxVelData():
    print "Max Velocity Data"
    print "Altitude:", u[numpy.argmax(u[:,1]),0]
    print "Velocity:", max(u[:,1])
    print "Time:",u[numpy.argmax(u[:,1]),2]
    

def getMinPosData():
    print "Min Position Data"
    print "Altitude:", u[numpy.argmin(u[:,0])-1,0]
    print "Velocity:", u[numpy.argmin(u[:,0])-1,1]
    print "Time:",u[numpy.argmin(u[:,0])-1,2]
    

getPosVel()
getMaxAltData()
getMaxVelData()
getMinPosData()
plotAll()
