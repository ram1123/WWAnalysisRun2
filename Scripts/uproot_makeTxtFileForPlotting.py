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
temp = []
for sample in ymload:
  temp.append(sample+".root")
  print "====================> ",sample
  nEvents=0
  nNegEvents=0
  for i,files in enumerate(Arrayfilepath):
    if files.find(sample) != -1:
    	print files
	#GetnEvents, nNegEvents = uproot.open(source+"/"+files)["otree"]["nEvents","nNegEvents"]
	print "#"*51,"\nDEBUG: 1"
	if not os.path.isfile(source+"/"+files): continue; 
	print "#"*51,"\nDEBUG: 2"
	if not (uproot.open(source+"/"+files).keys()):
		#exit(0)
		print "skip file: ",files
		print "#"*51,"\nDEBUG: 3"
	else:	
	  print "#"*51,"\nDEBUG: 4"
	  otree = uproot.open(source+"/"+files)["otree"]
	  InputArrays = otree.arrays(["nEvents","nNegEvents"])
	  #if GetnEvents and nNegEvents:
	  if otree:
    	  	temp.append(files)
    	  	if files.find('Single') != -1:
    	  		nEvents=1
    	  		nNegEvents=0
    	  	else:
    	  		if ifhaddOnly != 1:
    	  			#nEvents += int(GetnEvents.numentries)
	  			nEvents += InputArrays["nEvents"][0]
    	  			nNegEvents+= InputArrays["nNegEvents"][0] 
    	  		else:
    	  			print "skip..."
    	  	print files,nEvents,nNegEvents
	  #GetnEvents.close()

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
	#sampleInfo.sort()
	pprint(sampleInfo)
	
	#OutPutFile = "DibosonBoostedElMuSamples13TeV_"+datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M')+".txt";
	OutPutFile = "uproot_DibosonBoostedElMuSamples13TeV_"+source.split("/")[-4]+"_"+source.split("/")[-3]+"_"+source.split("/")[-2]+".txt";
	outScript = open(OutPutFile,"w");
	outScript.write("# name           file_location  xspb/lumipb  otherscale nMCevents       nMCNegEvents    colorcode       stackit\n")
	for i,files in enumerate(sampleInfo):
		#for sample in crossSections.keys():
		for sample in ymload:
			#print ymload[sections]
			#print sample
			#print i,files
			#print files,sample
			if files.find(sample) != -1:
				#print "===> ",files,sample
				print ymload[sample]["name"],"\t",files,"\t",ymload[sample]["CrossSection"],"\t1\t",int(ListnEvents[i]),"\t",int(ListnNegEvents[i]),"\t",ymload[sample]["ColorCode"],"\t",ymload[sample]["StackIt"]
				outScript.write(ymload[sample]["name"]+"\t"+files+"\t"+str(ymload[sample]["CrossSection"])+"\t1\t"+str(ListnEvents[i])+"\t"+str(ListnNegEvents[i])+"\t"+str(ymload[sample]["ColorCode"])+"\t"+str(ymload[sample]["StackIt"]))

stop = timeit.default_timer()
print('Time: ', stop - start)  
