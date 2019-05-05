#!/usr/bin/python
import os
import re
import sys
import yaml
import timeit
import datetime
import uproot 
# use of uproot package. This improves the running time by more than factor of 10
from colorpicker import * 

start = timeit.default_timer() # start timer

###################################################
#
#   User Inputs
#
###################################################
ifhaddOnly = 0    # Only hadd not text file then put ifhaddOnly=1 else ifhaddOnly=0

StoreArea = "/store/user/rasharma/SecondStep/Zmumu_2017_BugIssue_9X/WWTree_2019_05_04_03h33/"

PlottingDirectoryPath = "/uscms_data/d3/rasharma/aQGC_analysis/PlottingMacros/CMSSW_9_0_1/src/PlottingCodes2017/ControlPlots/SampleFiles/"

searchString = ['# name', 'data', 'Data', 'WV_EWK', 'aQGC', 'Diboson', 'VV', 'W\+jets', 'Z\+jets', 'top', 'QCD'] # this 

StoreAreaHadd = StoreArea+'HaddedFiles_Test/'
source = "/eos/uscms"+StoreArea
OutPutDir = source + 'HaddedFiles_Test'

#def CheckDirectory(path, directoryName):
#    """ Check if directory exists then delete and create a new one.
#    """
#    if os.path.isdir(path+directoryName):
#      print(CRED+"Directory "+source+directoryName+' found. Delete it...'+CEND)
#      #eos root://cmseos.fnal.gov rm
#      print "#"*51
#      print "# list all file in main directory: Just to check if directory HaddedFiles_Test does not exists"
#      os.system('ls /eos/uscms'+path)
#      print "#"*51
#      os.system('eos root://cmseos.fnal.gov/ rm -r '+source+'HaddedFiles_Test')
#      os.system('eos root://cmseos.fnal.gov/ mkdir '+source+'HaddedFiles_Test/')
#    else:
#      os.system('eos root://cmseos.fnal.gov/ mkdir '+source+'HaddedFiles_Test/')

# Check if the Directory exists, if yes then delete and then create else create new one.
if os.path.isdir(OutPutDir):
  print(CRED+"Directory "+source+'HaddedFiles_Test'+' found. Delete it...'+CEND)
  #eos root://cmseos.fnal.gov rm
  print "#"*51
  print "# list all file in main directory: Just to check if directory HaddedFiles_Test does not exists"
  os.system('ls /eos/uscms'+StoreArea)
  print "#"*51
  os.system('eos root://cmseos.fnal.gov/ rm -r '+source+'HaddedFiles_Test')
  os.system('eos root://cmseos.fnal.gov/ mkdir '+source+'HaddedFiles_Test/')
else:
  os.system('eos root://cmseos.fnal.gov/ mkdir '+source+'HaddedFiles_Test/')


Arrayfilepath = []
for root, dirs, filenames in os.walk(source):
	for f in filenames:
		filepath = f
		if filepath.endswith(".root"):
			Arrayfilepath.append(filepath)

print(Arrayfilepath)
print("====================================================")

List = []
ListnEvents = []
ListnNegEvents = []
print "#"*51

file = open("DataMCInfo.yml","r")
ymload = yaml.load(file)
file.close()

for sample in ymload:
  temp = []
  temp.append(sample+".root")
  print(CBOLD+"====================> "+sample+CEND)
  nEvents=0
  nNegEvents=0
  for i,files in enumerate(Arrayfilepath):
    if files.find(sample) != -1:
      if os.path.isfile(source+"/"+files): 
        if (uproot.open(source+"/"+files).keys()) == []:
        	print (CGREEN+"\nskip file: "+files+"\n"+CEND)
        else:	
	  #print (source+"/"+files)
          otree = uproot.open(source+"/"+files)["otree"]
          InputArrays = otree.arrays(["nEvents","nNegEvents"])
	  if len(InputArrays["nEvents"]):
            if files.find('Single') != -1:
              nEvents=1
              nNegEvents=0
            else:
              if ifhaddOnly != 1:
                nEvents += InputArrays["nEvents"][0]
                nNegEvents+= InputArrays["nNegEvents"][0] 
              else:
                print "skip..."
            temp.append(files)
            print files,nEvents,nNegEvents
      else:
        print(CRED+"\nFile Not Found: "+files+"\n"+CEND)

  if len(temp)>1:
  	List.append(temp)
  	ListnEvents.append(nEvents)
  	ListnNegEvents.append(nNegEvents)

print "#"*51
print(List)
print "#"*51


sampleInfo = []
for samples in List:
   temp=""
   temp="hadd -f "
   for i in range(0,len(samples)):
   	if i == 0:
   		temp+=source+'HaddedFiles_Test/'+samples[i]+' '
   	else:
   		temp+=source+"/"+samples[i]+" "
   print temp
   os.system(temp)
   sampleInfo.append(samples[0])
   print ""


if ifhaddOnly != 1:
  print "=============	MAKE SUMMARY	================"	
  print(sampleInfo)
  
  OutPutFile = "DibosonBoostedElMuSamples13TeV_"+source.split("/")[-4]+"_"+source.split("/")[-3]+"_"+source.split("/")[-2]+".txt";
  outScript = open(OutPutFile,"w");
  outScript.write("# name           file_location  xspb/lumipb  otherscale nMCevents       nMCNegEvents    colorcode       stackit\n")
  for i,files in enumerate(sampleInfo):
    for sample in ymload:
      if files.find(sample) != -1:
      	print ymload[sample]["name"],"\t",files,"\t",ymload[sample]["CrossSection"],"\t1\t",int(ListnEvents[i]),"\t",int(ListnNegEvents[i]),"\t",ymload[sample]["ColorCode"],"\t",ymload[sample]["StackIt"]
      	outScript.write(ymload[sample]["name"]+"\t"+StoreAreaHadd+files+"\t"+str(ymload[sample]["CrossSection"])+"\t1\t"+str(ListnEvents[i])+"\t"+str(ListnNegEvents[i])+"\t"+str(ymload[sample]["ColorCode"])+"\t"+str(ymload[sample]["StackIt"])+"\n")
  
  ################################################
  #
  #	Sort the file manually
  #
  ################################################
  outScript.close()

  SortedOutPutFile = "Sorted_"+OutPutFile
  outScript2 = open(SortedOutPutFile,'w')
  
  for strings in searchString:
    infile = open(OutPutFile,"r")
    for lines in infile:
      if re.match(str(strings),lines): 
        #print lines,
	outScript2.write(lines)
  outScript2.close()
  
################################################

# copy the created text file to the respective directory
os.system('xrdcp -f '+OutPutFile+'  root://cmseos.fnal.gov/'+ StoreAreaHadd)
os.system('xrdcp -f '+SortedOutPutFile+'  root://cmseos.fnal.gov/'+ StoreAreaHadd)
os.system('cp '+OutPutFile+' '+ SortedOutPutFile + ' ' +PlottingDirectoryPath)
# Stop timer
stop = timeit.default_timer()
# Print total time to run in minutes
print 'Run Time: ',round((stop - start)/60.0,2),'min'  
