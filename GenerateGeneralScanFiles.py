from GeneralCoordinates import *
import sys

maxRunTime = 1000000

if __name__ == '__main__':
    
    nsteps = int(sys.argv[1])
    setupType = sys.argv[2].upper()
    magnetList = [str(c).upper() for c in sys.argv[3].strip('[]').split(',')]
    bendName = sys.argv[4]
    momentum = float(sys.argv[5])
    entriesPerSplit = float(sys.argv[6])
    
    scanName = '/eos/user/c/camussol/SWAN_projects/OpticsStudies/Thesis/Coordinates/' + setupType
    
    
    if setupType in ['FDF', 'DFD', 'DDF', 'FFD', 'FDD', 'DFF']:
        
        nCoords = GenerateCoordinatesTriplets(scanName, nsteps, magnetList, momentum, entriesPerSplit)
                        
    elif setupType in ['FD', 'DF']:
                                
        nCoords = GenerateCoordinatesDoublets(scanName, nsteps, magnetList, momentum, entriesPerSplit)
    
    else:
        print('ERROR: Unknown setup type')
    
    
    
    
    if nCoords // entriesPerSplit != 0 or nCoords < entriesPerSplit:
        nFiles = int(nCoords/entriesPerSplit) + 1
        
    else:
        nFiles = int(nCoords/entriesPerSplit)
    
    
    lowerBound = 0
    upperBound = entriesPerSplit - 1
    
    subCommandTemplateName = 'SubmitTemplates/submissionTemplate.sh'
    subFileTemplateName = 'SubmitTemplates/submissionTemplateSub.sub'
    runCommandTemplateName = 'SubmitTemplates/runCommandTemplate.sh'
    
    subCommFileName = 'SubmitCommands/submission' + setupType
    subFileName = 'SubmitCommands/submission' + setupType
    runCommandName = 'SubmitCommands/runCommand' + setupType
        
    for i in range(nFiles):

        subCommandTemplate =  open(subCommandTemplateName, 'rt') 
        subFileTemplate = open(subFileTemplateName, 'rt') 
        runCommandTemplate =  open(runCommandTemplateName, 'rt') 
            
        with open(subCommFileName + str(i)+ '.sh', 'wt' ) as subFile:
            
            for line in subCommandTemplate :
                text = line %{'subFileName' : subFileName  + str(i)}
                subFile.write(text)

                
        with open(subFileName + str(i) + '.sub', 'wt') as subCommFile:
            for line in subFileTemplate :
                text = line %{'runCommandName' : runCommandName + str(i), 'splitID' : i, 'MaxRunTime' : maxRunTime, 'OutName' :  setupType + '_' + str(i), 'LogName' : setupType + '_' + str(i)}
                subCommFile.write(text)
           
                    
        with open(runCommandName + str(i) + '.sh', 'wt') as runFile:
            for line in runCommandTemplate :
                text = line %{'setupType' : setupType, 'bendName' : bendName, 'parameterFile' : scanName + str(i) + '.txt', 'scanStart' : lowerBound, 'scanEnd' : upperBound}
                runFile.write(text)
            
        lowerBound += entriesPerSplit 
        upperBound += entriesPerSplit
        
        if upperBound > nCoords:
            upperBound = nCoords
            
        subCommandTemplate.close()
        subFileTemplate.close()
        runCommandTemplate.close()
           
        
        
    
    
        
    
        
    
