#!/usr/bin/python
import uproot
import yaml
import sys
import datetime
import os
from pprint import pprint
import os
import ROOT as ROOT
import timeit

CRED = '\033[91m'
CGREEN  = '\33[32m'
CBOLD     = '\33[1m'
CEND = '\033[0m'

start = timeit.default_timer()

file = open("DataMCInfo.yml","r")
ymload = yaml.load(file)
file.close()

#print ymload

#source = "/eos/uscms/store/user/rasharma/SecondStep/WWTree_MuonPtScale_L1PreFire_2019_02_25_00h57/"
source = "/eos/uscms/store/user/rasharma/SecondStep/Run_2017/Frameworkupdate/WWTree_2019_04_04_14h15/"
os.system('xrdfs root://cmseos.fnal.gov/ mkdir '+source+'HaddedFiles/')

ifhaddOnly = 0

Arrayfilepath = []
for root, dirs, filenames in os.walk(source):
	for f in filenames:
		#filepath = root + os.sep + f
		filepath = f
		if filepath.endswith(".root"):
			Arrayfilepath.append(filepath)

pprint(Arrayfilepath)
pprint("====================================================")

#import difflib
ROOT.gROOT.SetBatch(True)

List = []
ListnEvents = []
ListnNegEvents = []
print "#"*51
for sample in ymload:
  temp = []
  temp.append(sample+".root")
  print(CBOLD+"====================> "+sample+CEND)
  nEvents=0
  nNegEvents=0
  for i,files in enumerate(Arrayfilepath):
    if files.find(sample) != -1:
      #print files
      if os.path.isfile(source+"/"+files): 
        #print "File to run: ",files
        if (uproot.open(source+"/"+files).keys()) == []:
        	print (CGREEN+"\nskip file: "+files+"\n"+CEND)
		#printline = "Skip file (No keys found): "+files
		#cprint(printline,'green')
        else:	
	  #print "Found keys: ",files
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
	    #GetnEvents.close()
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
   		temp+=source+'HaddedFiles/'+samples[i]+' '
   	else:
   		temp+=source+"/"+samples[i]+" "
   print temp
   os.system(temp)
   sampleInfo.append(samples[0])
   print ""

print "Ram Krishna Sharma"

if ifhaddOnly != 1:
	print "=============	MAKE SUMMARY	================"	
	print(sampleInfo)
	
	OutPutFile = "uproot_DibosonBoostedElMuSamples13TeV_"+source.split("/")[-4]+"_"+source.split("/")[-3]+"_"+source.split("/")[-2]+".txt";
	outScript = open(OutPutFile,"w");
	outScript.write("# name           file_location  xspb/lumipb  otherscale nMCevents       nMCNegEvents    colorcode       stackit\n")
	for i,files in enumerate(sampleInfo):
	  for sample in ymload:
	    if files.find(sample) != -1:
	    	print ymload[sample]["name"],"\t",files,"\t",ymload[sample]["CrossSection"],"\t1\t",int(ListnEvents[i]),"\t",int(ListnNegEvents[i]),"\t",ymload[sample]["ColorCode"],"\t",ymload[sample]["StackIt"]
	    	outScript.write(ymload[sample]["name"]+"\t"+files+"\t"+str(ymload[sample]["CrossSection"])+"\t1\t"+str(ListnEvents[i])+"\t"+str(ListnNegEvents[i])+"\t"+str(ymload[sample]["ColorCode"])+"\t"+str(ymload[sample]["StackIt"])+"\n")

stop = timeit.default_timer()
print('Time: ', stop - start)  
