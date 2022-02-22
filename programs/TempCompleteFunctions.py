# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 12:33:13 2021

@author: crazy
"""

import numpy as np
import scipy as scp
from scipy import signal
import matplotlib.colors as colors

def cenDiff(x,y): #4th order central difference stencil method
    '''
    x: first set of data
    y: second set of data y(x)
    
    return centered difference dydx
    '''
    if np.size(x) != np.size(y): #Check for simmilar size between x and y
        print("x and y not the same size!")
        return
    
    nx = np.size(x) #Get the number of x values
    dydx = np.zeros(nx) #Create an array of 0's with a simmilar size to x as stored in nx
    h = 0.001
    
    dydx[0] = (((-3*y[0]) + (4*y[1]) + (-1*y[2])) / (2*h)) / (((-3*x[0]) + (4*x[1]) + (-1*x[2])) / (2*h))
    dydx[1] = (((-3*y[0+1]) + (4*y[1+1]) + (-1*y[2+1])) / (2*h)) / (((-3*x[0+1]) + (4*x[1+1]) + (-1*x[2+1])) / (2*h))
    
    for i in range(2,nx-2):
        dydx[i] = (((y[i-2]) - (8*y[i-1]) + (8*y[i+1]) - (y[i+2])) / (12*h)) / (((x[i-2]) - (8*x[i-1]) + (8*x[i+1]) - (x[i+2])) / (12*h))
    
    dydx[-1] = (((3*y[nx-1]) + (-4*y[nx-2]) + (y[nx-3])) / (2*h)) / (((3*x[nx-1]) + (-4*x[nx-2]) + (x[nx-3])) / (2*h))
    dydx[-2] = (((3*y[nx-1-1]) + (-4*y[nx-2-1]) + (y[nx-3-1])) / (2*h)) / (((3*x[nx-1-1]) + (-4*x[nx-2-1]) + (x[nx-3-1])) / (2*h))
    
    return dydx #Return the dydx list

def plusGeneral(theta,phi,t,Qxx,Qyy,Qxy,Qxz,Qyz,Qzz,CenDiff=False):
    '''
    Get h+ in terms of arbitrary angle [M.A.P 4 June'19]
    If True, .dat files contain first time derivative of Q, else
    .dat files contain second time derivative of Q
    
    theta: theta information
    phi: phi angle information
    t: time information
    Q--: Quadrupole information (xx,yy,xy,xz,yz,zz)
    cenDiff: True (take centered difference), False (don't take centered difference)
    
    return h (gravitational wave information)
    
    Ken-ichi. Oohara, Takashi. Nakamura, Masaru. Shibata, Chapter 3. A Way to 3D Numerical Relativity, Progress of Theoretical Physics Supplement, Volume 128, March 1997, Pages 183–249, https://doi.org/10.1143/PTPS.128.183
    '''
    if (CenDiff == True): #If CenDiff == True Qddot arrays have yet to be differentiated
        QddotXX = cenDiff(t,Qxx)
        QddotYY = cenDiff(t,Qyy)
        QddotXY = cenDiff(t,Qxy)
        QddotXZ = cenDiff(t,Qxz)
        QddotYZ = cenDiff(t,Qyz)
        QddotZZ = cenDiff(t,Qzz)
    else:
        QddotXX = Qxx
        QddotYY = Qyy
        QddotXY = Qxy
        QddotXZ = Qxz
        QddotYZ = Qyz
        QddotZZ = Qzz
        
    '''Check these constants'''  
    Gc4 = 6.67e-8/3e10/3e10/3e10/3e10 #Constant
    distance = 10e3 * 3.086e18 # [10 kpc]
    
#     h = np.outer(QddotXX, np.cos(theta)*np.cos(theta)*np.cos(phi)*np.cos(phi)-np.sin(phi)*np.sin(phi)) + \
#         np.outer(QddotYY, np.cos(theta)*np.cos(theta)*np.sin(phi)*np.sin(phi)-np.cos(phi)*np.cos(phi)) - \
#         np.outer(QddotXY, np.sin(2*phi)*np.sin(theta)*np.sin(theta)) + \
#         np.outer(QddotZZ, np.sin(theta)*np.sin(theta)) - \
#         np.outer(QddotXZ, np.sin(2*theta)*np.cos(phi)) - \
#         np.outer(QddotYZ, np.sin(2*theta)*np.sin(phi))
    
    h = np.outer(QddotXX, np.cos(theta)*np.cos(theta)*np.cos(phi)*np.cos(phi)-np.sin(phi)*np.sin(phi)) + \
        np.outer(QddotYY, np.cos(theta)*np.cos(theta)*np.sin(phi)*np.sin(phi)-np.cos(phi)*np.cos(phi)) + \
        np.outer(QddotXY, np.sin(2*phi)*np.cos(theta)*np.cos(theta)+np.sin(2*phi)) + \
        np.outer(QddotZZ, np.sin(theta)*np.sin(theta)) - \
        np.outer(QddotXZ, np.sin(2*theta)*np.cos(phi)) - \
        np.outer(QddotYZ, np.sin(2*theta)*np.sin(phi))
    
    h *= 1*Gc4/distance
    
    h = np.reshape(h,(np.shape(QddotYZ)[0],np.shape(theta[0])[0],np.shape(theta[0])[0]))
    
    return(h)

def crossGeneral(theta,phi,t,Qxx,Qyy,Qxy,Qxz,Qyz,CenDiff=False):
    '''
    Get hx in terms of arbitrary angle [M.A.P 4 June'19]
    If True, .dat files contain first time derivative of Q, else
    .dat files contain second time derivative of Q
    
    theta: theta information
    phi: phi angle information
    t: time information
    Q--: Quadrupole information (xx,yy,xy,xz,yz)
    cenDiff: True (take centered difference), False (don't take centered difference)
    
    return h (gravitational wave information)
    
    Ken-ichi. Oohara, Takashi. Nakamura, Masaru. Shibata, Chapter 3. A Way to 3D Numerical Relativity, Progress of Theoretical Physics Supplement, Volume 128, March 1997, Pages 183–249, https://doi.org/10.1143/PTPS.128.183
    '''
   
    if (CenDiff == True):
        QddotXX = cenDiff(t,Qxx)
        QddotYY = cenDiff(t,Qyy)
        QddotXY = cenDiff(t,Qxy)
        QddotXZ = cenDiff(t,Qxz)
        QddotYZ = cenDiff(t,Qyz)
    else:
        QddotXX = Qxx
        QddotYY = Qyy
        QddotXY = Qxy
        QddotXZ = Qxz
        QddotYZ = Qyz

    Gc4 = 6.67e-8/3e10/3e10/3e10/3e10
    distance = 10e3 * 3.086e18 # [10 kpc]
  
    h = np.outer(QddotYY-QddotXX,np.cos(theta)*np.sin(phi)*np.cos(phi)) + \
        np.outer(QddotXY,np.cos(theta)*np.cos(2*phi)) + \
        np.outer(QddotXZ,np.sin(theta)*np.sin(phi)) - \
        np.outer(QddotYZ,np.sin(theta)*np.cos(phi))

    h *= 1*Gc4/distance
    
    h = np.reshape(h,(np.shape(QddotYZ)[0],np.shape(theta[0])[0],np.shape(theta[0])[0]))
    return(h)

def normDiffGeneral(theta,phi,t,Qxx,Qyy,Qxy,Qxz,Qyz,Qzz,CenDiff=False):
    '''
    theta: theta information
    phi: phi angle information
    t: time information
    Q--: Quadrupole information (xx,yy,xy,xz,yz,zz)
    cenDiff: True (take centered difference), False (don't take centered difference)
    
    return h (gravitational wave information)
    '''
    hPlusData = plusGeneral(theta,phi,t,Qxx,Qyy,Qxy,Qxz,Qyz,Qzz,CenDiff=CenDiff)
    hCrossData = crossGeneral(theta,phi,t,Qxx,Qyy,Qxy,Qxz,Qyz,CenDiff=CenDiff)
    
    hNormDiff = (hPlusData - hCrossData) / np.sqrt(hPlusData**2 + hCrossData**2)    

    return(hNormDiff)

def strain3dSurfMaxValue(dataObj, polarization, angles = [0,np.pi/2,0,2*np.pi]):
        '''
        Get the x,y,z values of strain ploted onto a 3d sphere ??
        
        dataObj: the object containing the data being acted upon
        polarization: setting for which polarization to use, (0 = cross, 1 = plus, 2 = norm diff)
        angles: angles to use in h calculations ([theta_0,theta,phi_0,phi] Domain: Theta:[0,pi] Phi:[0,2pi])
        '''
        #theta - altitudinal angle [radians]
        #phi   - azimuthal angle  [radians]
        thetaPhi = np.mgrid[angles[0]:angles[1]:dataObj.nVerts, angles[2]:angles[3]:dataObj.nVerts] #Define angles to use when calculating h
        
        analysisTime = dataObj.rawTime[dataObj.timeSample]
        
        QddotXX = cenDiff(analysisTime,dataObj.rawData[0][dataObj.timeSample])
        QddotYY = cenDiff(analysisTime,dataObj.rawData[1][dataObj.timeSample])
        QddotXY = cenDiff(analysisTime,dataObj.rawData[2][dataObj.timeSample])
        QddotXZ = cenDiff(analysisTime,dataObj.rawData[3][dataObj.timeSample])
        QddotYZ = cenDiff(analysisTime,dataObj.rawData[4][dataObj.timeSample])
        QddotZZ = cenDiff(analysisTime,dataObj.rawData[5][dataObj.timeSample])
        
        if polarization == 0: 
            hVals = crossGeneral(thetaPhi[0],thetaPhi[1],analysisTime,QddotXX,QddotYY,QddotXY,QddotXZ,QddotYZ,CenDiff=False)

        elif polarization == 1:
            hVals = plusGeneral(thetaPhi[0],thetaPhi[1],analysisTime,QddotXX,QddotYY,QddotXY,QddotXZ,QddotYZ,QddotZZ,CenDiff=False)

        elif polarization == 2:
            hVals = normDiffGeneral(thetaPhi[0],thetaPhi[1],analysisTime,QddotXX,QddotYY,QddotXY,QddotXZ,QddotYZ,QddotZZ,CenDiff=False)

        else:
            print('Correct the value given for polarization (0 = cross, 1 = plus, 2 = norm diff)!!')
            return()
        
        '''
        thetaPhi[0] -> theta
        thetaPhi[1] -> phi
        '''
        xSph = hVals * np.sin(thetaPhi[0]) * np.cos(thetaPhi[1])
        ySph = hVals * np.sin(thetaPhi[0]) * np.sin(thetaPhi[1])
        zSph = hVals * np.cos(thetaPhi[0])

        return(xSph,ySph,zSph,hVals)

def loadDatQ(path, fName, columns):
    '''
    path: file path to data files, not including file name 
    fName: file name
    columns: columns to load data from
    
    return data as 1d numpy arrays
    '''
    fileName  = path+"/"+fName
    
    t,Qxx,Qyy,Qxy,Qxz,Qyz,Qzz = 0, 0, 0, 0, 0, 0, 0
    
    t,Qxx,Qyy,Qxy,Qxz,Qyz,Qzz = np.loadtxt(fileName, usecols = (columns), unpack=True)
    
    return(t,Qxx,Qyy,Qxy,Qxz,Qyz,Qzz)

def interpAngles(angle, time, lenVector):
    '''
    Interpolate Angles for frequency sampling
    Angle: Array of angles
    Time: Array containing data time
    lenVector: Frequency range
    
    return interpolated angles
    '''
    angleReshape = np.reshape(angle, (len(angle)))
    
    fNew = scp.interpolate.interp1d(time,angleReshape) #Get interpolated angles
    
    times = np.linspace(time[0], time[-1], int(lenVector)) #Get times that will be used to sample fNew
    
    return(fNew(times))

def noNan(dataArray, conVal = 0):
    '''
    Convert any nan values in an array to 0
    
    dataArray: array to check through
    conVal: value to convert nan values too
    
    return array with nan elements replaced with 0 (by default)
    '''
    if np.sum(dataArray)==0:
        return(dataArray)
    
    dataArray[np.isnan(dataArray)] = 0
    
    return(dataArray)

def genSpectOutputs(N, timeData, angle1 = 0, angle2 = 0):
    '''
    Generate spectrogram information for the dipole direction of GW
    
    N: Frequency range uppler limit ??
    timeData: original time data [s]
    Angle1: Direction data
    Angle2: Direction data (optional)
    
    return (f1, t1, Sxx1, t2, f2, Sxx2)
    '''
    #Probably a better whay to accomplish the following 2 lines
    timeArray = np.linspace(timeData[0], timeData[-1], int(N)) #Generate an evenly spaced array of times of length N
    freqSample = N/(timeArray[-1]) #Generate frequency sampling value

    f1, t1, Sxx1 = 0, 0, 0
    if np.sum(angle1) != 0:
        angle1 = noNan(angle1) #Convert any nan values in the array angle1 to 0

        angle1Spec = interpAngles(angle1, timeData, N) #Time dependent angles information (interpolated)
        f1, t1, Sxx1 = signal.spectrogram(angle1Spec, freqSample) #Generate Spectrogram info

    f2, t2, Sxx2 = 0, 0, 0
    if np.sum(angle2) != 0:
        angle2 = noNan(angle2) #Convert any nan values in the array angle2 to 0
        
        angle2Spec = interpAngles(angle2, timeData, N) #Time dependent angles information (interpolated)
        f2, t2, Sxx2 = signal.spectrogram(angle2Spec, freqSample) #Generate Spectrogram info
    
    return([f1, t1, Sxx1, f2, t2, Sxx2])

def maxStrainDipoleDirection(dataObj, polarization = 0, angles = [np.pi/2,np.pi/2,0,2*np.pi], norm = False):
    '''
    
    '''
    iterStart = 0 #Iteration start
    iterEnd = len(dataObj.rawTime) #Iteration end

    iterRange = np.arange(iterStart,iterEnd,dataObj.timeStep) #Create iter array
    
    sV_X, sV_Y, sV_Z, sV_R = strain3dSurfMaxValue(dataObj, polarization = polarization, angles = angles) #sV_# -> strain value array (X,Y,Z,R)
    
    sV_R_Pos = np.array([]) #Initialize empty array

    sV_X_Flat = sV_X.reshape((len(sV_X),len(sV_X[0])**2)) #Flatten cord array to 2d array 
    sV_Y_Flat = sV_Y.reshape((len(sV_Y),len(sV_Y[0])**2)) #Flatten cord array to 2d array 
    sV_Z_Flat = sV_Z.reshape((len(sV_Z),len(sV_Z[0])**2)) #Flatten cord array to 2d array 
    sV_R_Flat = sV_R.reshape((len(sV_R),len(sV_R[0])**2)) #Flatten radius array to 2d array 

    sV_R_Pos = np.argmax(sV_R_Flat,axis=1) #Get max R value at each time

    sV_X_MV = np.take_along_axis(sV_X_Flat, np.expand_dims(sV_R_Pos, axis=1),axis=1) #Get values associated with max R values at each time
    sV_Y_MV = np.take_along_axis(sV_Y_Flat, np.expand_dims(sV_R_Pos, axis=1),axis=1) #Get values associated with max R values at each time
    sV_Z_MV = np.take_along_axis(sV_Z_Flat, np.expand_dims(sV_R_Pos, axis=1),axis=1) #Get values associated with max R values at each time

#     sV_XYZ_MV = hstack((sV_X_MV,sV_Y_MV,sV_Z_MV)) #Stack cord arrays 
    
    if norm == True: 
        rPos = (sV_X_MV**2 + sV_Y_MV**2 + sV_Z_MV**2)**(1/2) #Get length of r
        
        sV_X_MV = sV_X_MV/rPos #Normalize x data
        sV_Y_MV = sV_Y_MV/rPos #Normalize y data
        sV_Z_MV = sV_Z_MV/rPos #Normalize z data

        
#     print('Sum y max values', np.sum(sV_Y_MV))
#     print('Sum x max values', np.sum(sV_X_MV))
#     print('Sum of y/x', np.sum(sV_Y_MV/sV_X_MV))
#     print('Sum of r', np.sum(rPos))
    
    theta = np.arccos(sV_Z_MV/(np.sqrt((sV_X_MV**2)+(sV_Y_MV**2)+(sV_Z_MV**2))))
    phi = np.arctan(sV_Y_MV/sV_X_MV)
#     phi = np.arccos(sV_X_MV/sV_R_Pos)
    
    #Find different trig function to calculate phi above
    
    return(theta, phi, dataObj.rawTime[iterRange])

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    '''
    https://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib
    
    Generate a new colormap that matches the colormap for the full data for a subset of data
    
    cmap: full colormap
    minval: minimum value of the new colormap
    maxval: maximum value of the new colormap
    '''
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap