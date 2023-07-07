import subprocess
import os 
import time

commandsPath = r'/afs/cern.ch/user/c/camussol/private/Optics/SingleFile/SubmitCommands'
batchFiles = [f for f in os.listdir(commandsPath) if  f.endswith('.sh')]
submitFiles = [f for f in batchFiles if "submission" in f]

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
