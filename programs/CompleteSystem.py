# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 23:48:40 2022

@author: crazy
"""
import numpy as np
import TempCompleteFunctions as tpc
import matplotlib.pyplot as plt
import matplotlib.colors as colors

class GWAnalysis:
    '''
    Tool for analizing gravitational waves
    '''
    
    def __init__(self, dataPath, dataFile, dataColumns, timeSample = 0, timeStep = 5, tBounce = 0, nVerts = 5):
        '''
        dataPath: relative path to the main repository of the data/program
        dataFile: file name
        dataColumns: list of columns numbers that reference the data in the file [t,Qxx,Qyy,Qxy,Qxz,Qyz,Qzz]
        timeSample: range of times to work through when analizing data, set as list [start time index, end time index]
        timeStep: time step to take when analizing data
        tBounce: post bounce time [s]
        nVerts: number of vertecies/points to use when interpreting angle information
        '''
        #Import neccesary packages to the correct namespace
        
        
        if isinstance(nVerts, complex) == False:
            raise ValueError('nVerts must be entered as a complex number (ex: 10j), please enter a new value!')
        
        if nVerts.imag < 10j.imag:
            raise ValueError('nVerts can not be less than 10, please enter a new value!')
        
        #Set attributes
        self.dataPath = dataPath
        self.dataFile = dataFile
        self.dataColumns = dataColumns
        self.timeSample = timeSample
        self.timeStep = timeStep
        self.tBounce = tBounce
        self.nVerts = nVerts
        
    
    
    def update_nVerts(self, nVerts):
        '''
        Update the value for nVerts
        
        nVerts: number of vertecies/points to use when interpreting angle information,
                must be complex number 10 or greater
        '''
        #Check if value is complex
        if isinstance(nVerts, complex) == False:
            raise ValueError('nVerts must be entered as a complex number (ex: 10j), please enter a new value!')
        
        #Check if value is greater than or equal to 10j
        if nVerts.imag < 10j.imag:
            raise ValueError('nVerts can not be less than 10, please enter a new value!')

        self.nVerts = nVerts
        
        print('The number of verticies or points has been updated to', self.nVerts)
    
    def update_timeSample(self, timeSample):
        '''
        Update the value for timeSample
        
        timeSample: range of times to work through when analizing data, set as list [start time, end time, step through]
        '''
        self.timeSample = timeSample
        
    def update_timeStep(self, timeStep):
        '''
        Update the value for timeStep
        
        timeStep: time step to take when analizing data
        '''
        self.timeStep = timeStep
        
        self.update_timeSample(np.arange(0,len(self.rawTime),self.timeStep)) #Update timeSample with new value for timeStep
        
        print('The timeStep has been updated to', self.timeStep)
        
    def load_PreProcessedAngularData(self, dataFolder = '', skipRows = 0, delimeter = '', dtype = float):
        '''
        Load and store data
        
        dataFolder: name of folder where data is stored
        '''
        #Check if the data file is stored in a folder or not, create the correct filePath
        if dataFolder != 0:
            filePath = self.dataPath + '\\' + dataFolder + '\\' + self.dataFile
        
        else:
            filePath = self.dataPath + '\\' + self.dataFile
        
        t,phi,theta,magnitude = np.loadtxt(filePath, usecols = (self.dataColumns), unpack=True, skiprows = skipRows, delimiter = delimeter, dtype = dtype) 
    
        #Load in the data associated with the file and the specified columns
        self.rawTime = t
        self.rawData = np.array([theta,phi,magnitude])
        
        #Update timeSample
        self.update_timeSample(np.arange(0,len(self.rawTime),self.timeStep))
    
    def load_QuadData(self, dataFolder = 0):
        '''
        Load and store data
        
        dataFolder: name of folder where data is stored
        '''
        #Check if the data file is stored in a folder or not, create the correct filePath
        if dataFolder != 0:
            filePath = self.dataPath + '\\' + dataFolder + '\\' + self.dataFile
        
        else:
            filePath = self.dataPath + '\\' + self.dataFile
        
        #Load in the data associated with the file and the specified columns
        t,Qxx,Qyy,Qxy,Qxz,Qyz,Qzz = np.loadtxt(filePath, usecols = (self.dataColumns), unpack=True)
    
        self.rawTime = t
        self.rawData = np.array([Qxx,Qyy,Qxy,Qxz,Qyz,Qzz])
        
        #Update timeSample
        self.update_timeSample(np.arange(0,len(self.rawTime),self.timeStep))
        
    def set_rawDataHeaders(self, dataHeaders):
        '''
        Save headers for the raw data
        
        dataHeaders: list of headers equal in length to rawData
        '''
        self.dataHeaders = dataHeaders
        
#     def strain3dSurfMaxValue(self, polarization, angles = [0,np.pi,0,2*np.pi]):
#         '''
#         Get the x,y,z values of strain ploted onto a 3d sphere ??
        
#         polarization: setting for which polarization to use, (0 = cross, 1 = plus, 2 = norm diff)
#         angles: angles to use in h calculations ([theta_0,theta,phi_0,phi] Domain: Theta:[0,pi] Phi:[0,2pi])
#         '''
#         #theta - altitudinal angle [radians]
#         #phi   - azimuthal angle  [radians]
#         thetaPhi = np.mgrid[angles[0]:angles[1]:self.nVerts, angles[2]:angles[3]:self.nVerts] #Define angles to use when calculating h
        
#         analysisTime = self.rawTime[self.timeSample]
        
#         QddotXX = tpc.cenDiff(analysisTime,self.rawData[0][self.timeSample])
#         QddotYY = tpc.cenDiff(analysisTime,self.rawData[1][self.timeSample])
#         QddotXY = tpc.cenDiff(analysisTime,self.rawData[2][self.timeSample])
#         QddotXZ = tpc.cenDiff(analysisTime,self.rawData[3][self.timeSample])
#         QddotYZ = tpc.cenDiff(analysisTime,self.rawData[4][self.timeSample])
#         QddotZZ = tpc.cenDiff(analysisTime,self.rawData[5][self.timeSample])
        
#         if polarization == 0: 
#             hVals = tpc.crossGeneral(thetaPhi[0],thetaPhi[1],analysisTime,QddotXX,QddotYY,QddotXY,QddotXZ,QddotYZ,CenDiff=False)

#         elif polarization == 1:
#             hVals = tpc.plusGeneral(thetaPhi[0],thetaPhi[1],analysisTime,QddotXX,QddotYY,QddotXY,QddotXZ,QddotYZ,QddotZZ,CenDiff=False)

#         elif polarization == 2:
#             hVals = tpc.normDiffGeneral(thetaPhi[0],thetaPhi[1],analysisTime,QddotXX,QddotYY,QddotXY,QddotXZ,QddotYZ,QddotZZ,CenDiff=False)

#         else:
#             print('Correct the value given for polarization (0 = cross, 1 = plus, 2 = norm diff)!!')
#             return()
        
#         '''
#         Get confirmation for these conversions
#         '''
#         xSph = hVals * np.sin(thetaPhi[0]) * np.cos(thetaPhi[1])
#         ySph = hVals * np.sin(thetaPhi[0]) * np.sin(thetaPhi[1])
#         zSph = hVals * np.cos(thetaPhi[0])

#         return(xSph,ySph,zSph,hVals)
           
    def GW2dAnalysis(self, polarization, angles = [0,np.pi,0,2*np.pi], graphSize = [10,10], xlim = [0,0], ylim = [0,0]):
        '''
        
        polarization: setting for which polarization to use, (0 = cross, 1 = plus, 2 = norm diff)
        angles: angles to use in h calculations ([theta_0,theta,phi_0,phi] Domain: Theta:[0,pi] Phi:[0,2pi])
        graphSize: size of the image produced, list [width, height]
        xlim: x-axis limits for the graph, list [left bound, right bound]
        ylim: y-axis limits for the graph, list [lower bound, upper bound]
        '''
        #Data acquisition
        tempH = tpc.strain3dSurfMaxValue(self, polarization, angles)[3] #Only grab hVals

        tempH = tempH[:,4,0]# simplification/generalization of tempH = tempH.reshape(len(self.timeSample))
        
        #Plotting
        fig, ax = plt.subplots(1)
        fig.set_size_inches(graphSize[0],graphSize[1])

#         print(self.rawTime[self.timeSample])
#         print(self.rawTime[self.timeSample]-self.tBounce)
        
        ax.plot(self.rawTime[self.timeSample], tempH) 
        ax.axvline(self.tBounce, c='black', label = 'Post Bounce Time: ' + str(self.tBounce), alpha = 0.5)

        #Check if the limits are set by the user or not
        if xlim != [0,0]:
            ax.set_xlim(xlim)
        else:
            ax.set_xlim([-0.05,self.rawTime[-1]-self.tBounce+0.05])

        if ylim != [0,0]:
            ax.set_ylim(ylim)
            
        ax.set_title('GW Strain Vs Post Bounce Time ('+  self.dataFile +')')
        ax.set_xlabel('Time (s)')
        ax.legend()
        
#         attachedDataFR = np.loadtxt('RawDataFiles/s40_fr_GW_short.txt')

#         ax.plot(attachedDataFR[:,0], attachedDataFR[:,1], color = 'red')

#         ax.set_xlabel('Time (s)')
        
        
    def genScatterPlot(self, cord_XYZ_Set, prefix, path, title, analysisTime, figureSize = (15,15), dpi = 300, rotationArray = [0], norm = False, dataSlice = (0,2), axesLimit = [0,0], scaleFactor = [1,1,1], markerSize = 60, rasterizeTF = False, alpha = 0.8): 
        '''
        Check if timeIt is needed here
            Check to see if its the same length as xSet
            
        cord_XYZ_Set: array of cordinates x,y,z with shape [3,n] where n is the number of samples
        prefix: generated image filename prefix
        path: path from working directory to where generated image files should be stored (end with /)
        title: title of the generated plot
        rotationArray: array of degrees to rotate through when generating the images for the animation
        analysisTime: 
        figureSize: size of ploted figure
        dpi: pixels per inch for ploted figure
        norm: conditional, whether to normalize the plotted data or not
        dataSlice: tuple indicating what times (between the values) from the data will be plotted
        axesLimit: set specific axis limits for x,y,z (equal)
        scaleFactor: factors by which to scale the x,y,z axes
        markerSize: size for scatter plot markers
        rasterizeTF: check whether to rasterize the saved image or not
        alpha: alpha level to plot data points at
        '''

        #Unpack cord values (using .flatten on sliced object)
        xSet, ySet, zSet = cord_XYZ_Set[:,0], cord_XYZ_Set[:,1], cord_XYZ_Set[:,2]

        if norm == True: 
            rPos = (xSet**2 + ySet**2 + zSet**2)**(1/2) #Get length of r

            xSet = xSet/rPos #Normalize x data
            ySet = ySet/rPos #Normalize y data
            zSet = zSet/rPos #Normalize z data

        count = 0
        newLen = len(xSet[dataSlice[0]:dataSlice[1]]) #Length of the data after slicing with dataSliceVals

        colorMapName = 'viridis'
        colorMap = plt.get_cmap(colorMapName) #cmap to use
        new_cmap = tpc.truncate_colormap(colorMap, dataSlice[0]/(len(xSet)),
                                    dataSlice[1]/(len(xSet)), newLen) #New cmap generated from dataSliceVals and colorMap
        
        for i in rotationArray: #Redraw plot and savefig
#             fig = plt.figure(figsize = figureSize, constrained_layout = True)
            fig = plt.figure(figsize = figureSize)
            axSV_MV = fig.add_subplot(projection='3d')

            xPlot = xSet[dataSlice[0]:dataSlice[1]]
            yPlot = ySet[dataSlice[0]:dataSlice[1]]
            zPlot = zSet[dataSlice[0]:dataSlice[1]]

            '''
            The following code snippet is used to mask parts of the data
            In this case only the positive elements along the x-axis
            '''
    #         maskXPositive = xPlot > 0

    #         xPlot = xPlot[maskXPositive]
    #         yPlot = yPlot[maskXPositive]
    #         zPlot = zPlot[maskXPositive]

    #         newLen = len(xPlot)

            sc = axSV_MV.scatter(xPlot, yPlot, zPlot,
                                 c = np.linspace(dataSlice[0], dataSlice[1], newLen), 
                                 cmap = new_cmap, s = markerSize, edgecolors = 'black', alpha = alpha
                                ,linewidth = 0.5)

            axSV_MV.view_init(elev=0, azim=i+45) #Change view parameters (spin along the azimuth)

            axSV_MV.set_xlabel('x') #Assign label
            axSV_MV.set_ylabel('y') #Assign label
            axSV_MV.set_zlabel('z') #Assign label
            
#             axSV_MV.set_xticks([])

#             axSV_MV.set_xlabel('')

#             axSV_MV.set_box_aspect((1,1,3))

            if axesLimit != [0,0]:
                axSV_MV.set_xlim(axesLimit)
                axSV_MV.set_ylim(axesLimit)
                axSV_MV.set_zlim(axesLimit)
        
            else:
                xLim = np.array(axSV_MV.get_xlim())

                axSV_MV.set_xlim(xLim * scaleFactor[0])
                axSV_MV.set_ylim(xLim * scaleFactor[1])
                axSV_MV.set_zlim(xLim * scaleFactor[2])
            
#             print(axSV_MV.get_w_lims())
            
            figNorm = colors.Normalize(self.rawTime[0],self.rawTime[-1]) # map time data to colorbar (https://bit.ly/3lnV0VR)
            cBar = fig.colorbar(plt.cm.ScalarMappable(norm = figNorm, cmap = colorMapName), ax = axSV_MV, shrink=0.25) #Add color bar to axes = axSV_MV
            cBar.set_label('Time [s]')

            
#             fig.subplots_adjust(left=0, right=0.95, bottom=0, top=1)
#             plt.tight_layout(pad = 0.1,rect=[0.24, 0.03, 0.8, 0.95])
            
#             print(len(xPlot))
        
            if rasterizeTF == True:
                fig.set_rasterized(True) #Rasterize Image
        
            i_str = str(count)
            suffix = i_str.rjust(4,'0')            
            '''
            Figure out how I want to deal with file extension
            '''
            
#             fig.savefig(path + prefix + suffix, dpi = dpi) #Format dealt with in global parameters
            fig.savefig(path + prefix + suffix, dpi = dpi, bbox_inches = 'tight') #Format dealt with in global parameters
            fig.clear()
            plt.close(fig)
            count += 1
        
    def maxStrainScatterPlot(self, prefix, path, outputFileName, title, figureSize = (15,15), dpi = 300, polarization = 0, angles = [0,np.pi,0,2*np.pi], norm = False, rotationArray = np.arange(0,90,1), dataSlice = (0,0), axesLimit = [0,0], scaleFactor = [1,1,1], rasterizeTF = False, markerSize = 60, returnPlotDataTF = False, plotTF = True):
        '''
        
        prefix: generated image filename prefix
        path: path from working directory to where generated image files should be stored (end with /)
        title: title of the generated plot, if 'none' is used no title will be generated
        figureSize: size of ploted figure
        dpi: pixels per inch for ploted figure
        polarization: setting for which polarization to use, (0 = cross, 1 = plus, 2 = norm diff)
        angles: angles to use in h calculations ([theta_0,theta,phi_0,phi] Domain: Theta:[0,pi] Phi:[0,2pi])
        norm: conditional, whether to normalize the plotted data or not
        rotationArray: array of degrees to rotate through when generating the images for the animation
        dataSlice: tuple indicating what times (between the values) from the data will be plotted
        axesLimit: set specific axis limits for x,y,z (equal)
        scaleFactor: factors by which to scale the x,y,z axes
        rasterizeTF: check whether to rasterize the saved image or not
        markerSize: size for scatter plot markers 
        returnPlotDataTF: conditonal on whether to return plotted data (xyz, time, dataSliceStart, dataSliceEnd)
        plotTF: conditional n whether to create and save a plot of the data
        '''               
        if dataSlice == (0,0):
            dataSlice = (self.rawTime[0],self.rawTime[-1])
        
        if title.lower() == 'none': #Check if title is defined or not
            title = ''
            
        else:
            title += ' (Data: '+ self.dataFile +')'

            if polarization == 0:
                title += ' (cross-polarized)'
            elif polarization == 1:
                title += ' (plus-polarized)'
            elif polarization == 2:
                title += ' (normalized difference)'
            else:
                title += ' (ERROR-polarization)'

            if norm == True:
                title += ' (normalized)'
            else:
                title += ' (un-normalized)'
            
        analysisTime = self.rawTime[self.timeSample] #Sampled time
        
        sV_X, sV_Y, sV_Z, sV_HR = tpc.strain3dSurfMaxValue(self, polarization, angles = [0,np.pi,0,2*np.pi]) #sV_# -> strain value array (X,Y,Z,R)
            
        #sV_HR_Pos = np.array([]) #Initialize empty array
        
        sV_X_Flat = sV_X.reshape((len(sV_X),len(sV_X[0])**2)) #Flatten cord array to 2d array 
        sV_Y_Flat = sV_Y.reshape((len(sV_Y),len(sV_Y[0])**2)) #Flatten cord array to 2d array 
        sV_Z_Flat = sV_Z.reshape((len(sV_Z),len(sV_Z[0])**2)) #Flatten cord array to 2d array 
        sV_HR_Flat = sV_HR.reshape((len(sV_HR),len(sV_HR[0])**2)) #Flatten radius array to 2d array 

        sV_HR_Pos = np.argmax(sV_HR_Flat,axis=1) #Get max R value at each time

        sV_X_MV = np.take_along_axis(sV_X_Flat, np.expand_dims(sV_HR_Pos, axis=1),axis=1) #Get values associated with max R values at each time
        sV_Y_MV = np.take_along_axis(sV_Y_Flat, np.expand_dims(sV_HR_Pos, axis=1),axis=1) #Get values associated with max R values at each time
        sV_Z_MV = np.take_along_axis(sV_Z_Flat, np.expand_dims(sV_HR_Pos, axis=1),axis=1) #Get values associated with max R values at each time

        sV_XYZ_MV = np.hstack((sV_X_MV,sV_Y_MV,sV_Z_MV)) #Stack cord arrays 

        dataSliceStart = int(np.where(analysisTime == analysisTime[(analysisTime >= dataSlice[0]) & (analysisTime <= dataSlice[1])][0])[0]) #Get slice parameters from dataSlice
        dataSliceEnd = int(np.where(analysisTime == analysisTime[(analysisTime >= dataSlice[0]) & (analysisTime <= dataSlice[1])][-1])[0]) #Get slice parameters from dataSlice
        
        if plotTF == True:
            self.genScatterPlot(sV_XYZ_MV, prefix, path, title, analysisTime, figureSize, dpi, rotationArray, norm, dataSlice = (dataSliceStart, dataSliceEnd), axesLimit = axesLimit, scaleFactor = scaleFactor, rasterizeTF = rasterizeTF, markerSize = markerSize) #Generate 3D scatter plot

        if title.lower() != '':
            print('ffmpeg compilation string: ', end="")
            print('ffmpeg -framerate 10 -i ' + prefix + '%04d.jpeg ' + outputFileName + 'Movie.mp4')
            
        if returnPlotDataTF == True:
            return(sV_XYZ_MV, analysisTime, dataSliceStart, dataSliceEnd)
        
        
    def dipoleDirecPlots(self, N = 0, polarization = 0, angles = [0,np.pi/2,0,2*np.pi], modifyAngles = '', freqLimitTheta = (0,3500), freqLimitPhi = (0,3500), figSize = (40,20), norm = False, phiTF = False, specTF = False, saveTf = False, fileName = ''):
        '''
        Generate spectrogram plots of data to analyze dipole direction of GW

        polarization: setting for which polarization to use, (0 = cross, 1 = plus, 2 = norm diff)
        angles: which angles to take the spectrogram over [theta_initial, theta_final, phi_initial, phi_final]
        modifyAngles: modify the angles passed into the function with the designated trig function, ('cos' = cosine, 'sin' = sine)
        figSize: size to create each image as
        norm: conditional for whether to normalize the data or not 
        N: number of sampled points
        phiTF: conditional for whether to generate plots with the relevant phi data
        specTF: conditional for whether to generate a spectrogram for the data
        saveTF: conditional for whether to save the image or not
        fileName: output file name used when saveTF = True, continues from self.path
        '''
        fntSize = 25 #Fontsize for all plots
        legendSize = 25 #Size of the legend (if I can get it working)    
        
        #Check modifyAngles for an input, apply modification on angles if so
        if modifyAngles.lower() == 'cos':
            angles = np.cos(angles)

        elif modifyAngles.lower() == 'sin':
            angles = np.sin(angles)
        
        elif modifyAngles.lower() != '': #Check if modifyAngles has a valid input
            raise print('ERROR: pass in a valid argument for modifyAngles')
            return
            
        theta, phi, time = tpc.maxStrainDipoleDirection(self, polarization = polarization, angles = angles, norm = norm) #Generate theta, phi, and time data to be used in the following calculations
  
        fig = plt.figure(figsize = figSize)
        axSV_MV_A1 = fig.add_subplot(411)
        
#         print('Theta angles summed', np.sum(theta))
        
        axSV_MV_A1.plot(self.rawTime[self.timeSample], theta, c = 'blue', label = r'$\theta$')
        axSV_MV_A1.tick_params(axis='both', which='major', labelsize=fntSize/2)
#         axSV_MV_A1.set_title(self.dataFile + ' : theta', fontsize = fntSize)
#         axSV_MV_A1.legend(prop = {'size':legendSize})
        axSV_MV_A1.grid()

        #Check for any modification on angles, update the figure title if so
        if modifyAngles.lower() == '': 
            axSV_MV_A1.set_title(self.dataFile + ' : theta', fontsize = fntSize)

        elif modifyAngles.lower() != '':
            axSV_MV_A1.set_title(self.dataFile + ' :' + modifyAngles + '(theta)', fontsize = fntSize)
            
        
        title = 'Fourier Transform of:' + self.dataFile
        
        #For poster plots
        self.Spec_theta = theta
        self.Spec_time = self.rawTime[self.timeSample]
            
        #Check if spectrograms are enabled, calculate and plot if so
        if specTF == True:
            f1, t1, Sxx1 = tpc.genSpectOutputs(N, self.rawTime[self.timeSample], angle1 = theta)[:3]

            axSV_MV_F1 = fig.add_subplot(412, sharex = axSV_MV_A1)
            
            #Check for any modification on angles, update the figure title if so
            if modifyAngles.lower() == '': 
                title = 'Spectrogram (' + self.dataFile + ' : theta)'
            
            elif modifyAngles.lower() != '':
                title = 'Spectrogram (' + self.dataFile + ' :' + modifyAngles + '(theta))'

            axSV_MV_F1.pcolormesh(t1, f1, Sxx1, shading='gouraud')
            axSV_MV_F1.grid()
            axSV_MV_F1.tick_params(axis='both', which='major', labelsize=fntSize/2)
            axSV_MV_F1.set_ylabel('Frequency [Hz]', fontsize = fntSize)
            axSV_MV_F1.set_xlabel('Time [sec]', fontsize = fntSize)
            axSV_MV_F1.set_title(title, fontsize = fntSize)
            axSV_MV_F1.set_ylim(freqLimitTheta[0],freqLimitTheta[1])
            
            #For poster plots
            self.Spec_t1 = t1
            self.Spec_f1 = f1
            self.Spec_Sxx1 = Sxx1

        #Check if phi is enabled, calculate if so
        if phiTF == True:

            if specTF == True:
                axSV_MV_A2 = fig.add_subplot(413)

            else:
                axSV_MV_A2 = fig.add_subplot(412)
            
            axSV_MV_A2.plot(self.rawTime[self.timeSample], phi, c = 'red', label = '$\phi$')
            axSV_MV_A2.tick_params(axis='both', which='major', labelsize=fntSize/2)
#             axSV_MV_A2.set_title(self.dataFile + ' : phi', fontsize = fntSize)
#             axSV_MV_A2.legend(prop = {'size':legendSize})
            axSV_MV_A2.grid()

            #Check for any modification on angles, update the figure title if so
            if modifyAngles.lower() == '': 
                axSV_MV_A2.set_title(self.dataFile + ' : phi', fontsize = fntSize)
            
            elif modifyAngles.lower() != '':
                axSV_MV_A2.set_title(self.dataFile + ' :' + modifyAngles + '(phi)', fontsize = fntSize)
            
            #Check if spectrograms are enabled, calculate and plot if so
            if specTF == True:
                f2, t2, Sxx2 = tpc.genSpectOutputs(N, self.rawTime[self.timeSample], angle2 = phi)[3:]

                axSV_MV_F2 = fig.add_subplot(414, sharex = axSV_MV_A2)

                #Check for any modification on angles, update the figure title if so
                if modifyAngles.lower() == '': 
                    title = 'Spectrogram (' + self.dataFile + ' : phi)'
            
                elif modifyAngles.lower() != '':
                    title = 'Spectrogram (' + self.dataFile + ' :' + modifyAngles + '(phi))'

                axSV_MV_F2.pcolormesh(t2, f2, Sxx2, shading='gouraud')
                axSV_MV_F2.grid()
                axSV_MV_F2.tick_params(axis='both', which='major', labelsize=fntSize/2)
                axSV_MV_F2.set_ylabel('Frequency [Hz]', fontsize = fntSize)
                axSV_MV_F2.set_xlabel('Time [sec]', fontsize = fntSize)
                axSV_MV_F2.set_title(title, fontsize = fntSize)
                axSV_MV_F2.set_ylim(freqLimitPhi[0],freqLimitPhi[1])
        
        fig.tight_layout()
        
        if saveTf == True:

            fig.savefig(self.path + fileName + '.jpeg')
            
    def mollweideSpec(self, dataType = 'Quad', passTheta = 0, passPhi = 0, origin = (0,0), polarization = 0, angles = [0,np.pi/2,0,2*np.pi], figSize = (40,20), norm = False, dataSlice = (0,0), produceRotAxisTf = False, saveTf = False, fileName = ''):
        '''
        
        origin: setting which defines the origin of the graph (orgin = (RA,DEC) offsets [radians])
        polarization: setting for which polarization to use, (0 = cross, 1 = plus, 2 = norm diff)
        angles: which angles to take the spectrogram over [theta_initial, theta_final, phi_initial, phi_final]
        figSize: size to create each image as
        norm: conditional for whether to normalize the data or not
        dataSlice: tuple indicating what times (between the values) from the data will be plotted
        produceRotAxisTf: plot the rotation axis of the data
        saveTF: conditional for whether to save the image or not
        fileName: output file name used when saveTF = True, continues from self.path
        '''
        
        if dataSlice == (0,0):
            dataSlice = (self.rawTime[0],self.rawTime[-1])
        
        fig = plt.figure(figsize = (15,15))
        axSV_MV = fig.add_subplot(111, projection='mollweide')

        if dataType.lower() == 'quad':
            
            if produceRotAxisTf == False: #Check if produceRotAxisTf is False
                theta, phi, time = tpc.maxStrainDipoleDirection(self, polarization = polarization, angles = angles, norm = norm) #Generate theta, phi, and time data to be used in the following calculations
            
            elif produceRotAxisTf == True:
                theta, phi, xRotAxis, yRotAxis, zRotAxis, time = tpc.maxStrainDipoleDirection(self, polarization = polarization, angles = angles, norm = norm, returnCart=True) #Generate theta, phi, x, y, z, and time data to be used in the following calculations
                
                #Convert RotAxis coordinates to spherical coordinates
                thetaRotAxis = np.arccos(zRotAxis/(np.sqrt((xRotAxis**2)+(yRotAxis**2)+(zRotAxis**2))))
                phiRotAxis = np.arctan2(yRotAxis,xRotAxis)
                
        elif dataType.lower() == 'angular':
            theta = self.rawData[0]
            phi = self.rawData[1]
            time = self.rawTime
        
#         ra = -phi #azimuthal
#         dec = theta
        ra = -(np.remainder(phi+3*np.pi-origin[0],2*np.pi)-np.pi) #azimuthal, and setup so the origin can be moved
        dec = np.remainder(theta+1*np.pi-origin[1],1*np.pi)-np.pi/2 #altitudinal, shifted by 90 deg, and setup so the origin can be moved
        
        analysisTime = time #Sampled time
        
        dataSliceStart = int(np.where(analysisTime == analysisTime[(analysisTime >= dataSlice[0]) & (analysisTime <= dataSlice[1])][0])[0]) #Get slice parameters from dataSlice
        dataSliceEnd = int(np.where(analysisTime == analysisTime[(analysisTime >= dataSlice[0]) & (analysisTime <= dataSlice[1])][-1])[0]) #Get slice parameters from dataSlice
        
        ra = ra[dataSliceStart:dataSliceEnd]
        dec = dec[dataSliceStart:dataSliceEnd]
        
        
            
        newLen = len(theta[dataSliceStart:dataSliceEnd]) #Length of the data after slicing with dataSliceVals
        
        colorMapName = 'viridis'
        colorMap = plt.get_cmap(colorMapName) #cmap to use
        new_cmap = tpc.truncate_colormap(colorMap, dataSliceStart/(len(theta)),
                                    dataSliceEnd/(len(theta)), newLen) #New cmap generated from dataSliceVals and colorMap
        
        axSV_MV.scatter(ra, dec, c = np.linspace(dataSliceStart, dataSliceEnd, newLen), 
                                 cmap = new_cmap, s = 60, edgecolors = 'black', alpha = 0.8) #generate the scatter plot over a mollweide projection
        
        #Check if produceRotAxisTf is True, if so convert the rotation axis coordinates to ra and dec, slice it, and plot it.
        if produceRotAxisTf == True: 
            raRotAxis = -(np.remainder(phiRotAxis+3*np.pi-origin[0],2*np.pi)-np.pi) #azimuthal, and setup so the origin can be moved
            decRotAxis = np.remainder(thetaRotAxis+1*np.pi-origin[1],1*np.pi)-np.pi/2 #altitudinal, shifted by 90 deg, and setup so the origin can be moved
            
            raRotAxis = raRotAxis[dataSliceStart:dataSliceEnd]
            decRotAxis = decRotAxis[dataSliceStart:dataSliceEnd]
            
            print(raRotAxis == ra)
            
            axSV_MV.scatter(raRotAxis, decRotAxis, c = np.linspace(dataSliceStart, dataSliceEnd, newLen), 
                                 cmap = new_cmap, s = 60, marker = 'x', edgecolors = 'black', alpha = 0.7) #generate the scatter plot over a mollweide projection
        
        axSV_MV.grid(True) #enable the plots grid
        
        figNorm = colors.Normalize(self.rawTime[0],self.rawTime[-1]) # map time data to colorbar (https://bit.ly/3lnV0VR)
        cBar = fig.colorbar(plt.cm.ScalarMappable(norm = figNorm, cmap = colorMapName), ax = axSV_MV, shrink=0.5) #Add color bar to axes = axSV_MV
        cBar.set_label('Time(s)', fontsize = 15)
        
        return(axSV_MV)