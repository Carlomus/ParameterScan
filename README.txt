########How the various files work########

On mechanics of code will only talk about doublet but triplet is identical 

GenerateGeneralScanFiles.py takes as input ====> 

-nsteps
-setup (one of FD, FFD, FDF etc..)
-Magnet list
-Bend Name
-Max momentum 
-entries per split

and returns====>

-Nsteps Coordinates files in eos
-Enough .sh files to scan the whole parameter space in files of max size of 'entries per split'.
-.sub files to submit
-.sh files to run the submits


After generatin files simply rum RunScans.py to submit everything. Logs are in Outputs/Logs 

RunScans.py does not take any arguments

RunScans.py returns ===>

- Outputs/Doublet_results_i_j.txt files with the results (spotsizes and acceptances) for rows i to j of parameter file
- Outputs/Doublet_parameters_i_j.txt files with the parameters used in rows i to j of parameter file (coordinates)
- Log files




The way it works under the hood:

GenerateDoubletScanFiles runs DoubletCoordinates which for loops parameters and writes them to a file
In DoubletCoordinates there are some parameteres that can be varied and may need to be:
-Min 'safety spacing' currently 4 cm
-Max distance between magnets - currently 5 m
-Min distance to final focus (in code later the min distance is this distance + the safety spacing + magnet spacing)
-Returns the total number of configurations

GenerateDoubletScanFiles then generatesall the necessary .sh and .sub files by filling in from the templates
The number of files and lines being considered is calculated in this file using an index and skipping rows to be read from coordinates file
Also used in this calculation is the total number of configurations retured by GenerateDoubletScanFiles

Finally, RunScans.py loops through all .sh files in the folder SubmitCommands, makes them excutable and then runs the ones containing the string 'submission'

-For spotsize calculations the distribution of particles must be assumed. It has been dodgily done so far using a gaussian weighting.
In this weighing the sigmas are not the same as real life AND the amplitude may need to be rethought as it might affect the std

NOTE!! Various things related to collimators (eg lenth and position) have been hard coded, be careful!!