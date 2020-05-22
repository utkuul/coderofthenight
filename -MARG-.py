# -*- coding: utf-8 -*-
"""
Created on Sat May 02 21:21:08 2020

@author: UTKU
"""
import serial
import time
import numpy as np
from pyquaternion import Quaternion
from visual.controls import *
import math

# System constants
deltat = 0.001; #sampling period
gyroMeasError = 2.7 # gyro measurement error in rad/s 
gyroMeasError_rad = gyroMeasError * 3.14159265358979/180; # in rad/s
gyroMeasDrift = 3.14159265358979 * (0.015 / 180); # gyroscope measurement error in rad/s/s (shown as ... deg/s/s)
beta = np.sqrt(3 / 4.) * gyroMeasError_rad; # compute beta
zeta = np.sqrt(3 / 4.) * gyroMeasDrift; # compute zeta
q1 = 1.0; q2 = 0.0; q3 = 0.0; q4 = 0.0; #initial
bx=1;bz=0; #reference direction of flux in earth frame
wbx = 0; wby = 0; wbz = 0 # estimate gyroscope biases error

scene.range = 5
scene.width = scene.height = 1920
scene.background = color.white
scene.forward = vector(0,0,-1)

xarrow=arrow(length=2, shaftwidth=.1,opacity =0,color=color.red,axis=vector(1,0,0))
yarrow=arrow(length=2, shaftwidth=.1,opacity =0,color=color.green,axis=vector(0,1,0))
zarrow=arrow(length=2, shaftwidth=.1,opacity =0, color=color.blue,axis=vector(0,0,1))

#xlabel = label(pos=(2.25,0,0), text='X', color=color.black,height =20,border =1,font ='sans')
#ylabel = label(pos=(0,2.25,0), text='Y',color=color.black,height =20,border =1,font ='sans')
#zlabel = label(pos=(0,0,2.25), text='Z',color=color.black,height =20,border =1,font ='sans')

#xxlabel = label(pos=(4.25,0,0), text='XX', color=color.black,height =20,border =1,font ='sans')
#yylabel = label(pos=(0,4.25,0), text='YY',color=color.black,height =20,border =1,font ='sans')
#zzlabel = label(pos=(0,0,4.25), text='ZZ',color=color.black,height =20,border =1,font ='sans')

yan=arrow(length=4,shaftwidth=.1,opacity =0,color=color.black,axis=vector(1,0,0))
ust=arrow(length=4,shaftwidth=.1,opacity =0,color=color.magenta,axis=vector(0,1,0))
yukari=arrow(length=4,shaftwidth=.1,opacity =0,color=color.orange,axis=vector(0,0,1))

label (pos=(0,2,0),  color=color.blue,text = 'MARG VISUALIZATION')
freezelabel = label(pos=(0,-2.,0), color=color.blue, text = "ROCK YOUR SENSOR")

Boxx = box(pos=vector(0,0,0), size=vector(1,.1,2),opacity =0.7, color=color.blue)

margData = serial.Serial('com3',115200)
time.sleep(1)

while(True):
  try:
     while(margData.inWaiting()==0):
            pass
     dataPack = margData.readline()
     dataPack = unicode(dataPack,'utf-8')
     splitPack = dataPack.split(',')
     ax = float(splitPack[0]); print 'Ax =',ax,'m/s**2 ,',
     ay = float(splitPack[1]); print 'Ay =',ay,'m/s**2 ,',
     az = float(splitPack[2]); print 'Az =',az,'m/s**2'
     wx = float(splitPack[3]); print 'Wx =',wx,'rad/s ,',
     wy = float(splitPack[4]); print 'Wy =',wy,'rad/s ,',
     wz = float(splitPack[5]); print 'Wz =',wz,'rad/s'
     mx = float(splitPack[6]); print 'Mx =',mx,'uT',
     my = float(splitPack[7]); print 'My =',my,'uT,',
     mz = float(splitPack[8]); print 'Mz =',mz, 'uT',"\n"

     #Madgwick Algorithm - for MARG
     Sacce = Quaternion(np.array([0,ax,ay,az])) #accelerometer measurements
     Sgyro = Quaternion(np.array([0,wx,wy,wz])) #gyroscope measurements
     Smag = Quaternion(np.array([0,mx,my,mz]))  #magnetometer measurements
     SEq = Quaternion(np.array([q1,q2,q3,q4]))  #estimation of orientation

     #Normalise the accelerometer measurement values
     norm = np.sqrt(ax**2 + ay**2 + az**2)
     ax /= norm; ay /= norm; az /= norm

     #Normalise the magnetometer measurement values
     norm = np.sqrt(mx**2 + my**2 + mz**2)
     mx /= norm; my /= norm; mz /= norm;

     #Auxiliary variables
     halfq1 = 0.5*q1; halfq2 =0.5*q2; halfq3 = 0.5*q3; halfq4 = 0.5*q4;
     twoq1 = 2*q1; twoq2 = 2*q2; twoq3 = 2*q3; twoq4 = 2*q4;
     twobx = 2*bx; twobz = 2*bz;
     twomx= 2*mx; twomz = 2*mz;

     twobxq1 = 2*bx*q1;twobxq2 = 2*bx*q2;twobxq3 = 2*bx*q3;twobxq4= 2*bx*q4;
     twobzq1 = 2*bz*q1;twobzq2 = 2*bz*q2;twobzq3 = 2*bz*q3;twobzq4 = 2*bz*q4;

     q1q2=q1*q2 ; q1q3 = q1*q3; q1q4 = q1*q4;
     q2q3=q2*q3 ; q2q4 = q2*q4; q3q4 = q3*q4;

     q2sq=q2*q2; q3sq =q3*q3; q4sq=q4*q4;

     #Compute the objective function
     f1 = twoq2 * q4 - twoq1 * q3 - ax
     f2 = twoq1 * q2 + twoq3 * q4 - ay
     f3 = 1 - (twoq2 * q2) - (twoq3 * q3) - az
     f4 = twobx * (0.5 - q3sq - q4sq) + twobz * (q2q4 - q1q3) - mx
     f5 = twobx * (q2q3 - q1q4) + twobz * (q1q2 + q3q4) - my
     f6 = twobx * (q1q3+ q2q4) + twobz * (0.5 - q2sq - q3sq) - mz
     f = np.array ([f1,f2,f3,f4,f5,f6])

     #Jacobian Matrix
     J = np.array ([[-twoq3,twoq4,-twoq1,twoq2],
                   [twoq2,twoq1,twoq4,twoq3],
                   [0,-2 * twoq2, -2 * twoq3,0],
                   [-twobzq3, twobzq4, -2*twobxq3-twobzq1, -2*twobxq4 + twobzq2],
                   [-twobxq4 + twobzq2, twobxq3 + twobzq1, twobxq2 + twobzq4, -twobxq1 + twobzq3],
                   [twobxq3, twobxq4 - 2*twobzq2, twobxq1 - 2*twobzq3, twobxq2 ]])


     #Compute and normalise the gradient = Jt.f
     qhatdot = J.T.dot(f)
     qhatdot /= np.linalg.norm(qhatdot)

     #print qhatdot[0],qhatdot[1] ,qhatdot[2],qhatdot[3]


     #compute angular estimated direction of the gyroscope error
     w_err_x = 2 * q1 * qhatdot[1] - 2*q2 * qhatdot[0] - 2*q3 * qhatdot[3] + 2*q4 * qhatdot[2]
     w_err_y = 2 * q1 * qhatdot[2] + 2*q2 * qhatdot[3] - 2*q3 * qhatdot[0] - 2*q4 * qhatdot[1]
     w_err_z = 2 * q1 * qhatdot[3] - 2*q2 * qhatdot[2] + 2*q3 * qhatdot[1] - 2*q4 * qhatdot[0]

     #print  w_err_x, w_err_y, w_err_z


     #compute and DESTROY the gyroscope biases
     wbx += w_err_x * deltat * zeta;
     wby += w_err_y * deltat * zeta;
     wbz += w_err_z * deltat * zeta;
     wx -= wbx;
     wy -= wby;
     wz -= wbz;

     #print  wbx,wby,wbz
     #print  wx,wy,wz


     #wb = np.array([wbx,wby,wbz])
     #w_err = np.array([w_err_x,w_err_y,w_err_z])
     #wb += w_err * deltat * zeta

     #Compute the quaternion derivative measured by the gyroscope
     qdotomega = 0.5*SEq*Sgyro

     #print qdotomega[0],qdotomega[1] ,qdotomega[2],qdotomega[3]


     #Compute,integrate and normalise the estimated quaternion derivative
     q1 += (qdotomega[0] - (beta*qhatdot[0]))*deltat
     q2 += (qdotomega[1] - (beta*qhatdot[1]))*deltat
     q3 += (qdotomega[2] - (beta*qhatdot[2]))*deltat
     q4 += (qdotomega[3] - (beta*qhatdot[3]))*deltat
     norm2 = np.sqrt(q1**2 + q2**2 + q3**2 + q4**2)
     q1 /= norm2;q2 /= norm2;q3 /= norm2;q4 /= norm2


     # compute flux in the earth frame
     hx = 2 * mx * (0.5 - q3sq - q4sq) + 2 * my * (q2*q3 - q1*q4) + 2 * mz * (q2*q4 + q1*q3)
     hy = 2 * mx * (q2*q3 + q1*q4) + 2 * my * (0.5 - q2sq - q4sq) + 2 * mz * (q3*q4 - q1*q2)
     hz = 2 * mx * (q2*q4 - q1*q3) + 2 * my * (q3*q4 + q1*q2) + 2 * mz * (0.5 - q2sq - q3sq)

     #normalise the flux vector to have only components in the x and z
     bx = np.sqrt((hx**2) + (hy**2));
     bz = hz;

     SEq = Quaternion(np.array([q1,q2,q3,q4]))
     #print SEq

     #euler
     pitch = -math.atan2(2*(q1*q2+q3*q4),1-2*(q2**2+q3**2))
     roll = math.asin(2*(q1*q3-q4*q2))
     yaw = -math.atan2(2*(q1*q4+q2*q3),1-2*(q3**2+q4**2)) - np.pi/2

     yyy=vector(np.sin(yaw)*np.cos(roll),np.cos(yaw)*np.cos(roll),np.sin(roll)) #yaw's are negated because of axis differences
     rate(1000)
     z=vector(0,0,1)
     xxx=np.cross(yyy,z)
     zzz=np.cross(xxx,yyy)
     zrot=zzz*np.cos(pitch) + np.cross(yyy,zzz)*np.sin(pitch) 
     yan.axis=yyy
     ust.axis=np.cross(yyy,zrot)
     yukari.axis=zrot

     Boxx.axis =yyy
     Boxx.up = zrot

     yan.length=4
     ust.length=4
     yukari.length=4
  except:
      pass
