#!/usr/bin/python

import sys
import os

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 2:
        print "please enter path of input directory as one of argumetns of command\n"
	print "python Move_FileToTopDirFromCrab.py  <Path_of_Top_Crab_Directory>\n"
        exit(0)


#source = sys.argv[1]

#source = "/eos/uscms/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples_New/" + sys.argv[1]+"/"
source = "/eos/uscms/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/" + sys.argv[1]+"/"
print "Path of input ROOT file : ",source
#os.system("eos root://cmseos.fnal.gov ls "+sys.argv[1])

search1='log'
search2='failed'
Arrayfilepath = []
ArrayfileName = []
for root, dirs, filenames in os.walk(source):   
        for f in filenames:
                filepath = root + os.sep + f
                if filepath.find(search2) == -1:        # Don't select file if it contains word failed in path
                        if filepath.find(search1) == -1:        # Don't select file if it contains wored log in path
                                if filepath.endswith(".root"):  # select only root files
                                        Arrayfilepath.append(filepath)
					ArrayfileName.append(f)
#print len(Arrayfilepath)                                       
#print Arrayfilepath[0]

os.system("rm temp_script.sh")

OutFileName = "temp_script.sh"

OutFile = open(OutFileName,'w')

for i in range(0,len(Arrayfilepath)):
        #print "====> Copying file:",i," => ",Arrayfilepath[i]
        #print("xrdcp -f "+ Arrayfilepath[i]+" "+ source)
	# eos root://cmseos.fnal.gov mv
        #print("eos root://cmseos.fnal.gov mv "+ Arrayfilepath[i]+" "+ source+"/"+ArrayfileName[i])
        #OutFile.write("eos root://cmseos.fnal.gov mv "+ Arrayfilepath[i]+" "+ source+"/"+ArrayfileName[i]+'\n')
        OutFile.write("xrdcp -f "+ Arrayfilepath[i]+" "+ source+"/"+ArrayfileName[i]+'\n')
        #os.system("eos root://cmseos.fnal.gov mv "+ Arrayfilepath[i]+" "+ source)
	#print('')

OutFile.write("\n\n\n")
OutFile.write('echo "-------------------------"\n')
OutFile.write('echo "ls '+source+' "\n')
OutFile.write('echo "eosrm -r '+source+'/Samples "\n')
OutFile.close()
print "Run command: source ",OutFileName
os.system('bash '+str(OutFileName))
#os.system('rm '+str(OutFile))
