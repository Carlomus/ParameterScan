import shutil


toRerun = {
    
    'FD' : [0],
    'DF' : [],
    'FFD': [],
    'DFF': [],
    'FDF': [],
    'DDF': [0, 1, 2, 4, 7, 8, 9, 18, 20, 21, 22, 24, 25, 28, 30, 31, 33, 34, 37, 41, 42, 46, 49, 53, 56, 57, 58, 64, 66, 68, 74, 76, 81, 82, 85, 86, 88, 97],
    'FDD': [0, 2, 3, 4, 5, 8, 9, 11, 12, 13, 15, 16, 17, 18, 19, 21, 23, 24, 25, 26, 28, 29, 30, 32, 33, 35, 36, 39, 41, 42, 43, 44, 45, 47, 48, 50, 53, 54, 55, 56, 58, 59, 64, 66, 68, 69, 71, 72, 75, 76, 78, 79, 85, 89, 90, 91, 92, 94, 95, 97, 98],
    'DFD': [8, 31, 35, 37, 41, 76],
    
}



for setup in list(toRerun.keys()):

    for entry in toRerun[setup]:
        original = r'/afs/cern.ch/user/c/camussol/private/Optics/SingleFile/SubmitCommands/submission%s.sh'%(setup + str(entry))
        target = r'/afs/cern.ch/user/c/camussol/private/Optics/SingleFile/ResubmitCommands/submission%s.sh'%(setup + str(entry))

        shutil.copyfile(original, target)



import subprocess
import os 
import time


commandsPath = r'/afs/cern.ch/user/c/camussol/private/Optics/SingleFile/ResubmitCommands'
batchFiles = [f for f in os.listdir(commandsPath) if  f.endswith('.sh')]
submitFiles = [f for f in batchFiles if ("submission" in f) ]

for file in batchFiles:
    filename = os.path.join(commandsPath, file)
    #print("chmod +x %s" %(filename))
    subprocess.check_output("chmod +x %s" %(filename), shell=True)
    
print("Made batch files")

for file in submitFiles:
        filename = os.path.join(commandsPath, file) 
        subprocess.Popen(["%s" %(filename)], shell=True)
        time.sleep(0.5)
        

#subprocess.Popen(["watch -n 5 condor_q -nobatch -dag"], shell=True)
