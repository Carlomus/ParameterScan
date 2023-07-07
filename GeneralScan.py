import numpy as np
from math import *
import pandas as pd
from scipy.sparse import coo_matrix, block_diag
import sys
import argparse


def CreateMagnet( gradient, diameter, coreLength, realLength):
    
    magnet = {
     'gradient'  : gradient,
     'aperture'  : diameter/2 * 10 ** -3,
     'coreLength': coreLength,
     'realLength': realLength, 
    }
    
    return(magnet)

AvailableMagnets = {
    'QWL' : CreateMagnet(19.5, 100, 2.948, 3.306),
    'QSL' : CreateMagnet(18.5, 100, 3.0, 3.29),
    'QTL' : CreateMagnet(24, 80, 2.990, 3.335),
    'QTS' : CreateMagnet(24, 80, 1.490, 1.835),
    'QFS' : CreateMagnet(18.9, 100, 0.8, 1.08),  #hard
    'QFL' : CreateMagnet(18.9, 100, 1.2, 1.481),  #hard
    'QDS' : CreateMagnet(19.2, 91, 0.82, 1.08),  #hard
    'QPL' : CreateMagnet(10.5, 200, 2, 2.46),
    'QPS' : CreateMagnet(10.5, 200, 1, 1.46),
    'QNL' : CreateMagnet(24, 80, 2.990, 3.314),
    'QSS' : CreateMagnet(20.6, 100, 1.20, 2.25), #forget about this
    'QCS' : CreateMagnet(49.5, 90, 0.99, 1.51), #forget about this
   
}


def CreateDipole(B, lBend, momentum, apertureX, apertureY):

    rho = 3.34 * momentum / B
    bendingMagnet = {'B' : B,
                    'lBend' : lBend,
                    'rho' : rho,
                    'theta' : 2 * np.arctan2(lBend,2 * rho),
                    'apertureX' : apertureX * 10 ** -3,
                    'apertureY' : apertureY * 10 ** -3,
                    }
    return(bendingMagnet)


BendingMagnets = {
    'StandardDipole' : CreateDipole(1.84, 2, 13, 150, 100),
    'MBP110' : CreateDipole(1.94, 2.0, 13, 150, 55),
    'MBP140' : CreateDipole(1.75, 2.0, 13, 150, 70),
    'MBP200' : CreateDipole(1.48, 2.0, 13, 150, 100)
}



#These create ___ functions return a dictionary with all the relevant information of the element

def Create_Focusing_Thick_Lense(k, l, aperturex, aperturey):
        
    Rmat = RMAT_Create_Focusing_Thick_Lense(k,l)
    element = {'Rmat' : Rmat,
              'Aperture': [aperturex, aperturey],
              'Type' : 'F Quad'}
    
    return(element)
    
    
def Create_Defocusing_Thick_Lense(k, l, aperturex, aperturey):
    
    Rmat = RMAT_Create_Defocusing_Thick_Lense(k,l)    
    element = {'Rmat' : Rmat,
              'Aperture': [aperturex, aperturey],
              'Type' : 'D Quad'}
    
    return(element)
    
def Create_Bend(rho, theta, aperturex, aperturey, alpha1 = None, alpha2 = None):

    Rmat = RMAT_Create_Bend(rho, theta, alpha1, alpha2)
    
    element = {'Rmat' : Rmat,
              'Aperture': [aperturex, aperturey],
              'Type' : 'Bend'}
    
    return(element)


def Create_Drift(length, aperturex, aperturey):
    
    Rmat = RMAT_Create_Drift(length)

    element = {'Rmat' : Rmat,
              'Aperture': [aperturex, aperturey],
              'Type' : 'Drift'}
    
    return(element)

def Create_Collimator(length, aperturex, aperturey):
    
    Rmat = RMAT_Create_Drift(length)

    element = {'Rmat' : Rmat,
              'Aperture': [aperturex, aperturey],
              'Type' : 'Collimator'}
    
    return(element)

#Create an R matrix in the two planes for a thick lense 
def RMAT_Create_Focusing_Thick_Lense(k, l): 
    k = abs(k)
    sqrtk = k ** 0.5
    
    R11 = np.cos(l * sqrtk )
    R12 = sqrtk**-1  * np.sin(l *   sqrtk )
    R21 = -sqrtk  * np.sin(l * sqrtk )
    R22 = np.cos(l * sqrtk )
    R33 = np.cosh(l * sqrtk )
    R34 = sqrtk**-1  * np.sinh(l *   sqrtk )
    R43 = sqrtk  * np.sinh(l *   sqrtk )
    R44 = np.cosh(l *  sqrtk )


    Rmat = np.matrix([[R11, R12, 0,   0,   0,   0],
                      [R21, R22, 0,   0,   0,   0],
                      [0,   0,   R33, R34, 0,   0],
                      [0,   0,   R43, R44, 0,   0],
                      [0,   0,   0,   0,   1,   0],
                      [0,   0,   0,   0,   0,   1]])
    
    return(Rmat)

#Defocusing is same code as above but swap two blocks when putting into R matrix
def RMAT_Create_Defocusing_Thick_Lense(k, l): 
    k = abs(k)
    sqrtk = k ** 0.5

    R33 = np.cos(l * sqrtk )
    R34 = sqrtk**-1  * np.sin(l *   sqrtk )
    R43 = -sqrtk  * np.sin(l * sqrtk )
    R44 = np.cos(l * sqrtk )
    R11 = np.cosh(l * sqrtk )
    R12 = sqrtk**-1  * np.sinh(l *   sqrtk )
    R21 = sqrtk  * np.sinh(l *   sqrtk )
    R22 = np.cosh(l *  sqrtk )


    Rmat = np.matrix([[R11, R12, 0,   0,   0,   0],
                      [R21, R22, 0,   0,   0,   0],
                      [0,   0,   R33, R34, 0,   0],
                      [0,   0,   R43, R44, 0,   0],
                      [0,   0,   0,   0,   1,   0],
                      [0,   0,   0,   0,   0,   1]])
    
    return(Rmat)




#Create an R matrix for a drift lense 
def RMAT_Create_Drift(length):

    Rmat = np.matrix([[1, length, 0,   0,      0,   0],
                      [0,   1,    0,   0,      0,   0],
                      [0,   0,    1,   length, 0,   0],
                      [0,   0,    0,   1,      0,   0],
                      [0,   0,    0,   0,      1,   0],
                      [0,   0,    0,   0,      0,   1]])
    
    return (Rmat)


#theta is the bend angle, alpha is the angle of beam wrt magnet normal. Defaults to theta/2
def RMAT_Create_Bend(rho, theta, alpha1 = None, alpha2 = None):
    
    if alpha1 == None or alpha2 == None:
        alpha1 = theta/2
        alpha2 = theta/2
        
    R11 = np.cos(theta-alpha1)/np.cos(alpha1)
    R12 = rho* np.sin(theta)
    R21 = -1 * (1- np.tan(alpha1)*np.tan(alpha2)) * np.sin(theta - alpha1 - alpha2) / (rho * np.cos(alpha1 + alpha2))
    R22 = np.cos(theta-alpha2)/np.cos(alpha2)
    
    R16 = rho*(1-np.cos(theta))
    R26 =  np.sin(theta) + (1 - np.cos(theta))* np.tan(alpha2)
    
    R33 = 1 - theta* np.tan(alpha1)
    R34 = theta * rho
    R43 = -(np.tan(alpha1) + np.tan(alpha2))/ rho
    R44 = 1 - theta * np.tan(alpha2)
    
    
    Rmat = np.matrix([[R11, R12, 0,   0,   0, R16],
                      [R21, R22, 0,   0,   0, R26],
                      [0,   0,   R33, R34, 0,   0],
                      [0,   0,   R43, R44, 0,   0],
                      [0,   0,   0,   0,   1,   0],
                      [0,   0,   0,   0,   0,   1]])
    
    
    return(Rmat)




#Calculate functions multiply everything out. DOublet and Triplet are self explenatory, hard coded to multiply out certain elements
# From list allows for more flexibility

def CalculateDoublet(mag1, mag2, drift1, drift2, drift3):
    I = np.eye(6)
    partialR = [I]
    apertureList = []
    elementTypes = []
    
    #Multiply out the elements such that you get all the R matrices after each element
    for element in [drift1, mag1, drift2, mag2, drift3]:
        
        partialR.append(element['Rmat'] * partialR[-1])
        apertureList.append(element['Aperture'])
        elementTypes.append(element['Type'])
    
    #Dont forget to drop the initial identity
    return(partialR[1:], apertureList, elementTypes)
        
def CalculateTriplet(mag1, mag2, mag3, drift1, drift2, drift3, drift4):
    I = np.eye(6)
    partialR = [I]
    apertureList = []
    elementTypes = []
    
    for element in [drift1, mag1, drift2, mag2, drift3, mag3, drift4]:
        partialR.append(element['Rmat'] * partialR[-1])
        apertureList.append(element['Aperture'])
        elementTypes.append(element['Type'])
        
    return(partialR[1:], apertureList, elementTypes)


#Doublet + 4 bend achro
def CalculateDoubletAcromat(quad1, quad2, bend1, bend2, bend3, bend4, drift1, drift2, drift3, drift4, drift5):
    I = np.eye(6)
    partialR = [I]
    apertureList = []
    elementTypes = []
    
    for element in [drift1, quad1, drift2, quad2, drift3, bend1, drift4, bend2, drift5, bend3, drift4, bend4]:
        partialR.append(element['Rmat'] * partialR[-1])
        apertureList.append(element['Aperture'])
        elementTypes.append(element['Type'])
        
    return(partialR[1:], apertureList, elementTypes)

def CalculateLineFromList(elementsList):
    I = np.asmatrix(np.eye(6))
    partialR = [I]
    apertureList = []
    elementTypes = []
    
    for element in elementsList:
        partialR.append(element['Rmat'] * partialR[-1])
        apertureList.append(element['Aperture'])
        elementTypes.append(element['Type'])
        
    return(partialR[1:], apertureList, elementTypes)


def ChooseSetup(setupType, variables, apertureBx, apertureBy, AvailableMagnets):
    
    if setupType in ['FDF', 'DFD', 'DDF', 'FFD', 'FDD', 'DFF']:
    
        magnetUsed1 = variables[0]
        magnetUsed2 = variables[1]
        magnetUsed3 = variables[2]
        k1 = float(variables[3])
        k2 = float(variables[4])
        k3 = float(variables[5])
        l1 = float(variables[6])
        l2 = float(variables[7])
        l3 = float(variables[8])
        l4 = float(variables[9])

        magnet1 = AvailableMagnets[magnetUsed1]
        magnet2 = AvailableMagnets[magnetUsed2]
        magnet3 = AvailableMagnets[magnetUsed3]
        lq1 = magnet1['coreLength']
        lq2 = magnet2['coreLength']
        lq3 = magnet3['coreLength']

        maxAperture1 = magnet1['aperture']
        maxAperture2 = magnet2['aperture']
        maxAperture3 = magnet3['aperture']

        if setupType[0] == 'F':
            q1 = Create_Focusing_Thick_Lense(k1, lq1, maxAperture1, maxAperture1)
            
        else:
            q1 = Create_Defocusing_Thick_Lense(k1, lq1, maxAperture1, maxAperture1)
            k1 = -k1
            
        if setupType[1] == 'F':
            q2 = Create_Focusing_Thick_Lense(k2, lq2, maxAperture2, maxAperture2)
            
        else:    
            q2 = Create_Defocusing_Thick_Lense(k2, lq2, maxAperture2, maxAperture2)
            k2 = -k2
            
        if setupType[2] == 'F':
            q3 = Create_Focusing_Thick_Lense(k3, lq3, maxAperture3, maxAperture3)
            
        else:
            q3 = Create_Defocusing_Thick_Lense(k3, lq3, maxAperture3, maxAperture3)
            k3 = -k3
        
        
        d1 = Create_Drift(l1, maxAperture1, maxAperture1)
        d2 = Create_Drift(l2, maxAperture2, maxAperture2)
        d3 = Create_Drift(l3, maxAperture3, maxAperture3)
        d4 = Create_Drift(l4, apertureBx, apertureBy)

        elementsList = [d1, q1, d2, q2, d3, q3, d4]
        
        paramsDict = {
            'magnetUsed1': magnetUsed1,
            'magnetUsed2': magnetUsed2,
            'magnetUsed3': magnetUsed3,
            'k1': k1,
            'k2': k2,
            'k3': k3,
            'l1': l1,
            'l2': l2,
            'l3': l3,
            'l4': l4,
            'lFirstPart': l1 + l2 + l3 + l4 + lq1 + lq2 + lq3
        }
        
        
        return(elementsList, paramsDict)
        
    elif setupType in ['FD', 'DF']:
        
        
        magnetUsed1 = variables[0]
        magnetUsed2 = variables[1]
        k1 = float(variables[2])
        k2 = float(variables[3])
        l1 = float(variables[4])
        l2 = float(variables[5])
        l3 = float(variables[6])


        magnet1 = AvailableMagnets[magnetUsed1]
        magnet2 = AvailableMagnets[magnetUsed2]
        lq1 = magnet1['coreLength']
        lq2 = magnet2['coreLength']
        maxAperture1 = magnet1['aperture']
        maxAperture2 = magnet2['aperture']

            
                
        if setupType[0] == 'F':
            q1 = Create_Focusing_Thick_Lense(k1, lq1, maxAperture1, maxAperture1)
            
        else:
            q1 = Create_Defocusing_Thick_Lense(k1, lq1, maxAperture1, maxAperture1)
            k1 = -k1
            
        if setupType[1] == 'F':
            q2 = Create_Focusing_Thick_Lense(k2, lq2, maxAperture2, maxAperture2)
            
        else:    
            q2 = Create_Defocusing_Thick_Lense(k2, lq2, maxAperture2, maxAperture2)
            k2 = -k2
            
        d1 = Create_Drift(l1, maxAperture1, maxAperture1)
        d2 = Create_Drift(l2, maxAperture2, maxAperture2)
        d3 = Create_Drift(l3, apertureBx, apertureBy)
              
        elementsList = [d1, q1, d2, q2, d3]
        
        paramsDict = {
            'magnetUsed1': magnetUsed1,
            'magnetUsed2': magnetUsed2,
            'k1': k1,
            'k2': k2,
            'l1': l1,
            'l2': l2,
            'l3': l3,
            'lFirstPart': l1 + l2 + l3 + lq1 + lq2
        }
        
        return(elementsList, paramsDict)
    
    else:
        elementsList = []
        print('ERROR: Unknown setup type')
        return(None)
    


def CreateSaveDataFrames(setupType):
    
    if setupType in ['FDF', 'DFD', 'DDF', 'FFD', 'FDD', 'DFF']:
        
        results = pd.DataFrame(columns = ['BLID', 'Area x', 'Area y', 'Product Area', 'Length x', 'Height xp', 'Length y', 'Height yp', 'Distance to waist x', 'Beamspot size x', 'Beamspot size y', 'Product beamspot', 'Magnification'])
        parameters = pd.DataFrame(columns = ['BLID', 'Magnet 1' , 'Magnet 2', 'Magnet 3', 'Bending Magnet', 'k1', 'k2', 'k3', 'd1', 'd2', 'd3', 'd4', 'd5', 'padLength', 'Collimator Length', 'Total length'])

        
    elif setupType in ['FD', 'DF']:
                
        results = pd.DataFrame(columns = ['BLID', 'Area x', 'Area y', 'Product Area', 'Length x', 'Height xp', 'Length y', 'Height yp', 'Distance to waist x', 'Beamspot size x', 'Beamspot size y', 'Product beamspot', 'Magnification'])
        parameters = pd.DataFrame(columns = ['BLID', 'Magnet 1' , 'Magnet 2', 'Bending Magnet', 'k1', 'k2', 'd1', 'd2', 'd3', 'd4', 'padLength', 'Collimator Length', 'Total length'])

    else:
        print('ERROR: Unknown setup type')
        
    return(results, parameters)







#This function returns all the parameters we are really interested in. Takes in full R matrix in x, y and list of apertures at each point. 
#Also sigmas of distributions.

def FindResults(beamline, apertureList, elementTypes):
    
    #Generate the points along the target we are interested in
    
    testerPoint  = np.linspace(-30, 30, 61) * 10 ** -3
    
    xpLimsMax = []
    xpLimsMin = []
    ypLimsMax = []
    ypLimsMin = []
    
    for point in testerPoint:
    
        xpMax, xpMin, ypMax, ypMin = CalculatePrimeLimits(beamline, point, apertureList)
         
        xpLimsMax.append(xpMax)
        xpLimsMin.append(xpMin)
        
        ypLimsMax.append(ypMax)
        ypLimsMin.append(ypMin)
    
    #Convert the lists of maxima to arrays
    xpLimsMax = np.asarray(xpLimsMax)
    xpLimsMin = np.asarray(xpLimsMin)
    ypLimsMax = np.asarray(ypLimsMax)
    ypLimsMin = np.asarray(ypLimsMin)

    
    #Return arrays of only the correctly accepted phase space
    acceptedx  = testerPoint[(xpLimsMax >= xpLimsMin)]
    acceptedxp = xpLimsMax[(xpLimsMax >= xpLimsMin)]
    acceptedy  = testerPoint[(ypLimsMax >= ypLimsMin)]
    acceptedyp = ypLimsMax[(ypLimsMax >= ypLimsMin)]
    
    #Get the length of accepted phase space at x' = 0
    #Get the height of accepted phase space at x = 0
    if len(acceptedx)>0:
    	xlength = acceptedx[-1] - acceptedx[0]
    else:
    	xlength = 0

    if len(acceptedxp) == 0:
	xpheight = 0 
    else:
	xpheight = 2 * acceptedxp[int(len(acceptedxp)/2)]

    if len(acceptedy  )>0:
    	ylength = acceptedy[-1] - acceptedy[0]
    else:
    	ylength = 0

    if len(acceptedyp) == 0:
	ypheight = 0
    else:
	ypheight = 2 * acceptedyp[int(len(acceptedyp)/2)]

    
    #Get the difference between max and min. This then enables to check whether xpMax> xpMin and so up to where particles are accepted.
    diffX = xpLimsMax - xpLimsMin
    diffY = ypLimsMax - ypLimsMin
    
    #Work out Phase space area
    areaX, areaY = FindArea(diffX, diffY)
    
    #Find where the collimator is 
    collPosition = elementTypes.index('Collimator')
    
    
    ###############################################################
    ###############################################################
    ################## HARD CODED COLLIMATOR LENGTH!!! KEEP IN MIND
    ###############################################################
    ###############################################################
    
    
    collLen  = 1.1 #################
    dpad = 0.7


    #Work out the size in x and y of the image at middle of 'ideal collimator' and distance from end of dpad
    focusDistance, bmsptx, bmspty = FindPhaseSpaceSizeAtCollimatorVariable(beamline,  collPosition, testerPoint, xpLimsMax, xpLimsMin, testerPoint, ypLimsMax, ypLimsMin, collLen, dpad)
    
    
    return(areaX, areaY, focusDistance, bmsptx, bmspty, xlength, xpheight, ylength, ypheight)

    
def CalculatePrimeLimits(RList, x, apertureList):
    
    #create lists of max and min values particles can have
    upperBoundX = []
    lowerBoundX = []
    upperBoundY = []
    lowerBoundY = []
    
    
    #Loop over each R matrix element and the aperture at that point. (Both defined at the end of the element)
    for R, aperture in zip(RList, apertureList):
        
        aperturex = aperture[0]
        aperturey = aperture[1]
        
        #Work out the angle a particle at a given x has to have to reach the aperture at either the top (+) or bottom (-) of the element
        xppos = (aperturex  - R[0,0] * x)/ R[0,1]
        xpneg = (-aperturex - R[0,0] * x)/ R[0,1]        
        yppos = (aperturey  - R[2,2] * x)/ R[2,3]
        ypneg = (-aperturey - R[2,2] * x)/ R[2,3]
        
        
        #If R12 is positive an increase in x' leads to an increase in the final x (ie Xf = x*R11 + x'*R12 so dXf =  R12*dx'), so x' calculated here will be the largest possible x' at that pos
        #(and vice versa for -, where x' will be the smallest possible angle)
        
        if R[0,1]>0:
            upperBoundX.append(xppos)
            lowerBoundX.append(xpneg)            
        #If R12 is negative however, an increase in x' will lead to a decrease of Xf, meaning that at +aperture the angle is actually the smallest acceptable angle and vice versa fot -aperture  
        elif R[0,1]<0:
            lowerBoundX.append(xppos)
            upperBoundX.append(xpneg)
        else: 
            continue
            
            
        if R[2,3]>0:
            upperBoundY.append(yppos)
            lowerBoundY.append(ypneg)            
        #If R12 is negative however, an increase in x' will lead to a decrease of Xf, meaning that at +aperture the angle is actually the smallest acceptable angle and vice versa fot -aperture  
        elif R[2,3]<0:
            lowerBoundY.append(yppos)
            upperBoundY.append(ypneg)
        else: 
            continue
    
    
    #Take the smallest of the largest angles as this will be what imposes the limit and the same by taking the max of the lower bounds
    xp1 = min(upperBoundX)
    xp2 = max(lowerBoundX)
    
    yp1 = min(upperBoundY)
    yp2 = max(lowerBoundY)
    
    #Returns maximum and minimum x'
    
    #NOTE: There is no guarantee that xpMax is > xpMin, and actually there will always be a cross over point. All this means is that the acceptance from that x is 0 -- no angles can make it between
    # +A and -A
    return(xp1, xp2, yp1, yp2)
    
    

def FindArea(diffX, diffY):
    
    #As arguments we have the differences between x'Max and x'Min, meaning that if x'Max > x'Min then entry is positive
    #Cut out the negative parts of the distributions and integrate the leftover, ie only the area between xpMax and min
    #Note the dx = 1e-3 is because the sampling of the points was every one mm - if this changes this will have to change too.
    
    posPartX = np.maximum(diffX, 0)
    areaX = np.trapz(posPartX, dx=1e-3)
    
    posPartY = np.maximum(diffY, 0)
    areaY = np.trapz(posPartY, dx=1e-3)
    
    return(areaX, areaY)


def FindPhaseSpaceSizeAtCollimatorVariable(Rmatrices, positionCollimator, testerPointX, xpMax, xpMin, testerPointY, ypMax, ypMin, collLen, dpad):
    
    xImageMid = 10000000
    safetySpacing = 0.34
    smaplingSpacing = 0.01
    
    
    R = Rmatrices[positionCollimator-2]
    Rx = R[0:2, 0:2]
    Ry = R[2:4, 2:4]
    
    
    #Cut the various xpMax and xpMin to only where xpMax is > xpMin, ie where acceptance is not 0
    xpMaxNew = xpMax[(xpMax >= xpMin)] 
    xpMinNew = xpMin[(xpMax >= xpMin)] 
    ypMaxNew = ypMax[(ypMax >= ypMin)]
    ypMinNew = ypMin[(ypMax >= ypMin)] 

    #Return the points along the targets where there are particles that can reach the end of the line
    pointXNew = testerPointX[(xpMax >= xpMin)] 
    pointYNew = testerPointY[(ypMax >= ypMin)]
    
    #Find where to place collimator such that image is smallest at the center of it
    
    for d in np.arange(safetySpacing + collLen/2, collLen/2 + 2 * dpad - safetySpacing, smaplingSpacing):
    
        absXLim, absYLim = CalculateImageSize(d, Rx, Ry, pointXNew, xpMaxNew, xpMinNew, pointYNew, ypMaxNew, ypMinNew)

        if absXLim <= xImageMid:
            xImageMid = absXLim
            optimalCollDistance = d - collLen/2 - dpad
            
        else: break
    
    #Look at the size of the image at start and end of collimator, see which one is limiting and pick that as the definition of resolution
    
    absXLimStart, absYLimStart = CalculateImageSize(d - collLen/2, Rx, Ry, pointXNew, xpMaxNew, xpMinNew, pointYNew, ypMaxNew, ypMinNew)
    absXLimEnd, absYLimEnd = CalculateImageSize(d + collLen/2, Rx, Ry, pointXNew, xpMaxNew, xpMinNew, pointYNew, ypMaxNew, ypMinNew)
    
    if absXLimStart > absXLimEnd:
        xImage = absXLimStart
        yImage = absYLimStart
        
    else:
        xImage = absXLimEnd
        yImage = absYLimEnd
    
    
    
    return(optimalCollDistance , xImage, yImage)



def CalculateImageSize(d, Rx, Ry, pointXNew, xpMaxNew, xpMinNew, pointYNew, ypMaxNew, ypMinNew):
    
    Rdriftx = np.matrix([[1,   d],
                        [0,   1]]) 
    
    Rdrifty = np.matrix([[1,   d],
                        [0,   1]])

    RxSpot = Rdriftx * Rx
    RySpot = Rdrifty * Ry

    R11 = RxSpot[0,0]
    R12 = RxSpot[0,1]
    R33 = RySpot[0,0]
    R34 = RySpot[0,1]
    
    XCollMax = R11 * pointXNew + R12 * xpMaxNew
    XCollMin = R11 * pointXNew + R12 * xpMinNew

    YCollMax = R33 * pointYNew + R34 * ypMaxNew
    YCollMin = R33 * pointYNew + R34 * ypMinNew

    absXLim = max([max(XCollMax), max(XCollMin)])
    absYLim = max([max(YCollMax), max(YCollMin)])
    
    return(absXLim, absYLim)
    


collLen = 1.1
aptCollx = 0.1
aptColly = 0.1
padLength = 0.7
lAchromat = 5


if __name__ == '__main__':

    setupType = sys.argv[1].upper()
    bendName = sys.argv[2]
    infileName = sys.argv[3]
    parameterRangeStart = int(sys.argv[4])
    parameterRangeEnd = int(sys.argv[5])
    

    results, parameters = CreateSaveDataFrames(setupType)
    i = parameterRangeStart     
    
    with open(infileName, 'r') as infile:

        for line in infile:
            
            variables = line.split()

            dipole = BendingMagnets[bendName]
            rho = dipole['rho'] 
            theta = dipole['theta']
            apertureBx = dipole['apertureX']
            apertureBy = dipole['apertureY']
            lbend = dipole['lBend']
                        
            elementsList, paramsDict = ChooseSetup(setupType, variables, apertureBx, apertureBy, AvailableMagnets)
            
            
            dAchromat = Create_Drift(lAchromat, apertureBx, apertureBy)
            coll = Create_Collimator(collLen, aptCollx, aptColly)                        
                        
            b1 = Create_Bend(rho, theta, apertureBx, apertureBy)
            b2 = Create_Bend(-rho, -theta, apertureBx, apertureBy)
            b3 = Create_Bend(-rho, -theta, apertureBx, apertureBy)
            b4 = Create_Bend(rho, theta, apertureBx, apertureBy)
              
            dcap = Create_Drift(1.016, 0.1, 0.1)      
                
            dpad = Create_Drift(padLength, apertureBx, apertureBy)
            elementsListAchromat = [ b1, dAchromat, b2, dpad, coll, dpad, b3, dAchromat, b4, dcap]
                        
                
            fullElementsList = elementsList + elementsListAchromat
                
            R, apts, elementTypes = CalculateLineFromList(fullElementsList)
            accAx, accAy, optimalCollDist, bmsptx, bmspty, xlength, xpheight, ylength, ypheight = FindResults(R, apts, elementTypes)
                        
            
            lenTot = paramsDict['lFirstPart'] + 2* lAchromat + 2 * padLength + collLen + 4 * lbend 
            
            #Save the results to the original dictionaries
            
            
            if setupType in ['FDF', 'DFD', 'DDF', 'FFD', 'FDD', 'DFF']:
        
                params = [str(i), paramsDict['magnetUsed1'], paramsDict['magnetUsed2'], paramsDict['magnetUsed3'], bendName, paramsDict['k1'], paramsDict['k2'], paramsDict['k3'], paramsDict['l1'], paramsDict['l2'], paramsDict['l3'], paramsDict['l4'], lAchromat, padLength, collLen, lenTot]
        
            elif setupType in ['FD', 'DF']:
                
                params = [str(i), paramsDict['magnetUsed1'], paramsDict['magnetUsed2'], bendName, paramsDict['k1'], paramsDict['k2'], paramsDict['l1'], paramsDict['l2'], paramsDict['l3'], lAchromat, padLength, collLen, lenTot]
            
            #append by writing to a new index (len of list is always +1 compared to largest index)
            parameters.loc[len(parameters.index)] = params
            res = [str(i), accAx, accAy, accAx * accAy, xlength, xpheight, ylength, ypheight, optimalCollDist, bmsptx, bmspty, bmsptx * bmsptx, 2 * bmsptx/xlength ]
            results.loc[len(results.index)] = res

            i += 1

        nameParams = '/eos/user/c/camussol/SWAN_projects/OpticsStudies/Thesis/ScanOutputs/Parameters_' + setupType + '_' + str(parameterRangeStart) + '_' + str(parameterRangeEnd)+ '.txt'
        nameResults = '/eos/user/c/camussol/SWAN_projects/OpticsStudies/Thesis/ScanOutputs/Results_' + setupType + '_' + str(parameterRangeStart) + '_' + str(parameterRangeEnd) + '.txt'

        parameters.to_csv(nameParams, float_format = '%.8f', index = False)
        results.to_csv(nameResults,  float_format = '%.8f', index = False)