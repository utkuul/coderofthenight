# -*- coding: utf-8 -*-
"""
Created on Fri May 01 17:17:59 2020

@author: UTKU
"""
import serial
import time
import numpy as np
from pyquaternion import Quaternion
from visual.controls import * 
import math

'''System Constants'''
deltat = 0.01 #sampling period
gyroMeasError = 2.183 # shown as 5 deg/s
gyroMeasError_rad = gyroMeasError * np.pi/180 # for rad/s
beta = np.sqrt(3/4.)*gyroMeasError_rad
q1 = 1.0; q2 = 0.0; q3 = 0.0; q4 = 0.0; #initial

scene.range = 5
scene.width = scene.height = 1000
scene.background = color.white
scene.forward = vector(0,0,-1)

xarrow=arrow(length=2, shaftwidth=.1, color=color.red,axis=vector(1,0,0))
yarrow=arrow(length=2, shaftwidth=.1, color=color.green,axis=vector(0,1,0))
zarrow=arrow(length=2, shaftwidth=.1, color=color.blue,axis=vector(0,0,1))

xlabel = label(pos=(2.25,0,0), text='X', color=color.black,height =20,border =1,font ='sans')
ylabel = label(pos=(0,2.25,0), text='Y',color=color.black,height =20,border =1,font ='sans')
zlabel = label(pos=(0,0,2.25), text='Z',color=color.black,height =20,border =1,font ='sans')

xxlabel = label(pos=(4.25,0,0), text='XX', color=color.black,height =20,border =1,font ='sans')
yylabel = label(pos=(0,4.25,0), text='YY',color=color.black,height =20,border =1,font ='sans')
zzlabel = label(pos=(0,0,4.25), text='ZZ',color=color.black,height =20,border =1,font ='sans')

yan=arrow(length=4,shaftwidth=.1,color=color.black,axis=vector(1,0,0))
ust=arrow(length=4,shaftwidth=.1,color=color.magenta,axis=vector(0,1,0))
dik=arrow(length=4,shaftwidth=.1,color=color.orange,axis=vector(0,0,1))

label (pos=(0,5.,0),  color=color.red,text = 'IMU VISUALIZATION')
freezelabel = label(pos=(0,-5.,0), color=color.red, text = "MOVE YOUR SENSOR")

Boxx = box(pos=vector(0,0,0), size=vector(1,.1,2), color=color.blue)

arduinoData = serial.Serial('com3',115200)
time.sleep(1)

while(True):
  try:  
     while(arduinoData.inWaiting()==0):
            pass
     dataPacket = arduinoData.readline()
     dataPacket = bytearray(dataPacket,'utf-8')
     splitPacket = dataPacket.split(',')
     ax = float(splitPacket[0]); print 'Ax =',ax,'m/s**2 ,',
     ay = float(splitPacket[1]); print 'Ay =',ay,'m/s**2 ,',
     az = float(splitPacket[2]); print 'Az =',az,'m/s**2'
     wx = float(splitPacket[3]); print 'Wx =',wx,'rad/s ,',
     wy = float(splitPacket[4]); print 'Wy =',wy,'rad/s ,',
     wz = float(splitPacket[5]); print 'Wz =',wz,'rad/s',"\n"
           
     #Madgwick Algorithm - for IMU
     Sacce = Quaternion(np.array([0,ax,ay,az])) #accelerometer measurements
     Sgyro = Quaternion(np.array([0,wx,wy,wz])) #gyroscope measurements
     SEq = Quaternion(np.array([q1,q2,q3,q4]))  #estimation of orientation
    
     #Normalise the accelerometer measurement values
     norm = np.sqrt(ax**2 + ay**2 + az**2)
     ax = ax/norm; ay = ay/norm; az = az/norm
    
     #Auxiliary variables 
     twoq1 = 2*q1; twoq2 = 2*q2; twoq3 = 2*q3; twoq4 = 2*q4;
    
     #Compute the objective function 
     f1 = twoq2 * q4 - twoq1 * q3 - ax
     f2 = twoq1 * q2 + twoq3 * q4 - ay
     f3 = 1 - (twoq2 * q2) - (twoq3 * q3) - az
     f = np.array ([f1, f2, f3]) # objective function

     #f = np.array ([
     #twoq2 * q4 - twoq1 * q3 - ax,
     #twoq1 * q2 + twoq3 * q4 - ay,
     #1 - (twoq2 * q2) - (twoq3 * q3) - az]) # objective function
    
     #Jacobian Matrix
     J = np.array ([[-twoq3,twoq4,-twoq1,twoq2],
                    [twoq2, twoq1, twoq4, twoq3],
                    [0,-2*twoq2, -2*twoq3,0]])
    
     #Compute and normalise the gradient = Jt.f
     qhatdot = J.T.dot(f) #8.basamakta yuvarliyor.
     qhatdot /= np.linalg.norm(qhatdot)
     
     #Compute the quaternion derivative measured by the gyroscope
     qdotomega = 0.5*SEq*Sgyro
    
     #Compute,integrate and normalise the estimated quaternion derivative 
     q1 += (qdotomega[0] - (beta*qhatdot[0]))*deltat
     q2 += (qdotomega[1] - (beta*qhatdot[1]))*deltat
     q3 += (qdotomega[2] - (beta*qhatdot[2]))*deltat
     q4 += (qdotomega[3] - (beta*qhatdot[3]))*deltat
     norm2 = np.sqrt(q1**2 + q2**2 + q3**2 + q4**2)
     q1 /= norm2;q2 /= norm2;q3 /= norm2;q4 /= norm2
     SEq = Quaternion(np.array([q1,q2,q3,q4]))
       
     #conversion of Quaternion to Euler
     roll= math.atan2(2*(q1*q2+q3*q4),1-2*(q2**2+q3**2))  
     pitch= -math.asin(2*(q1*q3-q4*q2))
     yaw= -math.atan2(2*(q1*q4+q2*q3),1-2*(q3**2+q4**2)) - np.pi/2

     #the rotation provided by basic trigonometry for pitch&yaw
     yyy=vector(np.sin(yaw)*np.cos(roll),np.cos(yaw)*np.cos(roll),np.sin(roll)) 
     rate(100)
     z=vector(0,0,1) 
     xxx=np.cross(yyy,z) 
     zzz=np.cross(xxx,yyy)

     #Rodrigues Rotation Formula 
     zrot=zzz*np.cos(pitch) + np.cross(yyy,zzz)*np.sin(pitch)    

     #imaginary axis
     yan.axis=yyy
     ust.axis=np.cross(yyy,zrot)
     dik.axis=zrot 
     Boxx.axis =yyy
     Boxx.up = zrot
        
     yan.length=4
     ust.length=4
     dik.length=4
  except:
     pass

