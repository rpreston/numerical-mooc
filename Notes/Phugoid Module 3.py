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
ve = 325
Cd = 0.15
mpo = 100
dt = 0.1
N = 800
z = numpy.zeros(N)
v = numpy.zeros(N)
t = numpy.zeros(N)
c = numpy.zeros(N)

#u = numpy.array([h0,v0])

def stepFunction(currentTime,startTime,endTime,highValue,lowValue):
    if currentTime>=startTime and currentTime<endTime :
        return highValue
    else:
        return lowValue


def getPropellantMass(currentTime):
    if currentTime<5 :
        returnValue = (mpo-currentTime*stepFunction(currentTime,0.,5.,20.,0.))
    else:
        returnValue = 0
    return returnValue

def getPosVel():
    currentTime = 0.0
    currentH = 0
    count = 1
    u = numpy.array([h0,v0])
    z[0] = h0
    while currentH>=0.0 :
        u = u + dt*numpy.array([u[1],(-g) + (((stepFunction(currentTime,0.,5.,20.,0.)*ve)-(1/2*rho*u[1]*numpy.absolute(u[1])*maxCrossSection*Cd)))/(totalRocketMass+getPropellantMass(currentTime))])
        
        z[count] = u[0]
        v[count] = u[1]
        t[count] = currentTime
        c[count] = count
        currentH = u[0]
        currentTime+=dt
        count += 1

def getMaxAltData():
    countNew = 0
    for item in z:
        if(item == max(z)):
            maxAltC = countNew
        countNew+=1
    print countNew
    print "Max Altitude Data"
    print "Altitude:", z[maxAltC]
    print "Velocity:",v[maxAltC]
    print "Time:",t[maxAltC]

def getMaxVelData():
    countNew = 0
    for item in v:
        if(item == max(v)):
            maxVelC = countNew
        countNew+=1
    print "Max Velocity Data"
    print "Altitude:", z[maxVelC]
    print "Velocity:", v[maxVelC]
    print "Time:",t[maxVelC]
    
getPosVel()
getMaxAltData()
getMaxVelData()