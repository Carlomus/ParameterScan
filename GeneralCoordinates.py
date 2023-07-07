import numpy as np
import itertools
import sys

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



#Additional spacing on top of the minimum distance between magnets
extraSpacing = 0.400

maxDistanceDrifts = 5

minimumSpacingBend = 0.335

distanceQPLMBPL = 0.451 + 0.230 + 0.335  ####connector + end of qpl + end of mbpl

#Use these in the various k range definitions if you know what you want to start with
minklow  =  0.05
minkhigh =  0.12


def GenerateCoordinatesTriplets(outfileName, nsteps, magnetList, momentum, entriesPerSplit):
       
    fileID = 0
    j = 0
    
    outfile = open(outfileName + str(fileID) + ".txt", "w")
    
    
    #Loop over the two magnets wanted
    #for magnetUsed1, magnetUsed2, magnetUsed3 in set(itertools.product(magnetList, magnetList, magnetList)): 
    
    #Change the for loop to get all combinations of magnets. With bottom loop only get  same type throughout
    for magnetType in magnetList:       
        
        magnetUsed1 = magnetType 
        magnetUsed2 = magnetType 
        magnetUsed3 = magnetType 

        #Returns the dict of magnet parameters, which are then called for calculations
        magnetUsed1 = magnetUsed1.upper()
        magnetUsed2 = magnetUsed2.upper()
        magnetUsed3 = magnetUsed3.upper()

        magnet1 = AvailableMagnets[magnetUsed1]
        magnet2 = AvailableMagnets[magnetUsed2]
        magnet3 = AvailableMagnets[magnetUsed3]
        lq1 = magnet1['coreLength']
        lq2 = magnet2['coreLength']
        lq3 = magnet3['coreLength']
            
        maxk1 = magnet1['gradient'] * 0.2998 / momentum
        maxAperture1 = magnet1['aperture']
        maxk2 = magnet2['gradient'] * 0.2998 / momentum
        maxAperture2 = magnet2['aperture']
        maxk3 = magnet3['gradient'] * 0.2998 / momentum
        maxAperture3 = magnet3['aperture']

        #Additional spacing due to core lengths (where B field is) being shorter than actual magnet lengths
        minimumSpacing1 =  (magnet1['realLength'] - magnet1['coreLength'] ) / 2
        minimumSpacing2 =  (magnet2['realLength'] - magnet2['coreLength'] ) / 2
        minimumSpacing3 =  (magnet3['realLength'] - magnet3['coreLength'] ) / 2


        #Lists of parameters which we will scan through
        length1Range = np.linspace(minimumSpacing1  + extraSpacing , maxDistanceDrifts, nsteps)
        length2Range = np.linspace(minimumSpacing1 + minimumSpacing2 + extraSpacing , maxDistanceDrifts, nsteps)
        length3Range = np.linspace(minimumSpacing2 + minimumSpacing3  + extraSpacing, maxDistanceDrifts, nsteps)
        #length4Range = np.linspace(minimumSpacing3 + minimumSpacingBend, maxDistanceDrifts + 3, nsteps)
        length4Range = [distanceQPLMBPL]  #Dirty hack but d3 is fixed to be as small as possible
           
        #K values should technically start at 0, but since not very finely scanned we put an arbitrary start
        k1Range = np.linspace(minklow, maxk1, nsteps)
        k2Range = np.linspace(minklow, maxk2, nsteps)
        k3Range = np.linspace(minklow, maxk3, nsteps)

        

        #Loop over all parameters
        for k1Iter in k1Range:
            for k2Iter in k2Range:
                for k3Iter in k3Range:
                    for d1Iter in length1Range:
                        for d2Iter in length2Range:
                            for d3Iter in length3Range:
                                for d4Iter in length4Range:
                                    
                                    if j < entriesPerSplit:
                                            
                                        outfile.write("%s %s %s %f %f %f %f %f %f %f \n" %(magnetUsed1, magnetUsed2, magnetUsed3, k1Iter, k2Iter, k3Iter, d1Iter, d2Iter, d3Iter, d4Iter))
                                        j += 1
                                            
                                    else:
                                        
                                        outfile.close()
                                        fileID += 1
                                        j = 1
                                        outfile = open(outfileName + str(fileID) + ".txt", "w")
                                        outfile.write("%s %s %s %f %f %f %f %f %f %f \n" %(magnetUsed1, magnetUsed2, magnetUsed3, k1Iter, k2Iter, k3Iter, d1Iter, d2Iter, d3Iter, d4Iter))

    if not outfile.close():                                
        outfile.close()
    return(len(magnetList) * len(length1Range) * len(length2Range ) * len(length3Range) * len(length4Range ) * len(k1Range) * len(k2Range) * len(k3Range))





def GenerateCoordinatesDoublets(outfileName, nsteps, magnetList, momentum, entriesPerSplit):
        
    fileID = 0
    j = 0
        
    outfile = open(outfileName + str(fileID) + ".txt", "w")
    
    #Loop over the two magnets wanted
    #for magnetUsed1, magnetUsed2 in set(itertools.product(magnetList, magnetList)):        

    #Change the for loop to get all combinations of magnets. With bottom loop only get  same type throughout
    for magnetType in magnetList:       
        
        magnetUsed1 = magnetType 
        magnetUsed2 = magnetType 

        #Returns the dict of magnet parameters, which are then called for calculations
        magnetUsed1 = magnetUsed1.upper()
        magnetUsed2 = magnetUsed2.upper()

        magnet1 = AvailableMagnets[magnetUsed1]
        magnet2 = AvailableMagnets[magnetUsed2]
        lq1 = magnet1['coreLength']
        lq2 = magnet2['coreLength']
            
        maxk1 = magnet1['gradient'] * 0.2998 / momentum
        maxAperture1 = magnet1['aperture']
        maxk2 = magnet2['gradient'] * 0.2998 / momentum
        maxAperture2 = magnet2['aperture']

        #Additional spacing due to core lengths (where B field is) being shorter than actual magnet lengths
        minimumSpacing1 =  (magnet1['realLength'] - magnet1['coreLength'] ) / 2
        minimumSpacing2 =  (magnet2['realLength'] - magnet2['coreLength'] ) / 2


        #Lists of parameters which we will scan through
        length1Range = np.linspace(minimumSpacing1, maxDistanceDrifts, nsteps)
        length2Range = np.linspace(minimumSpacing1 + minimumSpacing2 + extraSpacing, maxDistanceDrifts, nsteps)
        length3Range = np.linspace(minimumSpacing2 + minimumSpacingBend  + extraSpacing, maxDistanceDrifts, nsteps)
        length3Range = [distanceQPLMBPL]  #Dirty hack but d3 is fixed to be as small as possible
        
        #K values should technically start at 0, but since not very finely scanned we put an arbitrary start
        k1Range = np.linspace(minkhigh, maxk1, nsteps)
        k2Range = np.linspace(minkhigh, maxk2, nsteps)
        
        #Loop over all parameters
        for k1Iter in k1Range:
            for k2Iter in k2Range:
                for d1Iter in length1Range:
                    for d2Iter in length2Range:
                        for d3Iter in length3Range:
                                
                            if j < entriesPerSplit:
                                
                                outfile.write("%s %s %f %f %f %f %f \n" %(magnetUsed1, magnetUsed2, k1Iter, k2Iter, d1Iter, d2Iter, d3Iter))
                                j +=1
                                    
                            else:
                                outfile.close()
                                fileID += 1
                                j = 1
                                outfile = open(outfileName + str(fileID) + ".txt", "w")
                                outfile.write("%s %s %f %f %f %f %f \n" %(magnetUsed1, magnetUsed2, k1Iter, k2Iter, d1Iter, d2Iter, d3Iter))
                                    
    if not outfile.close():                                
        outfile.close()                                
    return(len(magnetList) * len(length1Range) * len(length2Range ) * len(length3Range ) * len(k1Range) * len(k2Range))

