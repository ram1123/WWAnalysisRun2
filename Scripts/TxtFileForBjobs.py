#!/usr/bin/python

import sys
import os
import ROOT as ROOT
import subprocess

ROOT.gROOT.SetBatch(True)

from ConfigParser import RawConfigParser
config = RawConfigParser()   
config.optionxform = str       # Last two lines are done because ConfigParser will not preserve case
config.read("DataMCInfo.ini")

Name   = dict([sample, str(nPro)] for sample, nPro in config.items('FullName'))
crossSections = dict([sample, float(xsec)] for sample, xsec in config.items('CrossSection'))

#from pprint import pprint
#print "cross sections:" 
#pprint(Name)
source = "/eos/uscms/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/"

ArrayDirPath = []
for root, dirs, filenames in os.walk(source):
        #print dirs
        for d in dirs:
                #print d
                ArrayDirPath.append(root+os.sep+d)

OutputArray = [ ]

for sample in Name.keys():
        #print sample,Name[sample]
        nEvents = 0
        nNegEvents = 0
        nTotalEvents = 0
        nTotalNegEvents = 0
        for directory in sorted(ArrayDirPath):
                if directory.find(sample) != -1:
                        #print "==> ",sample,"\t:\t",directory
                        output = '1'
                        if directory.find('Single') != -1:
				print "This is data... Do nothing..."
                        else:
                                OutputArray.append("( "+str(crossSections[sample])+"\t"+directory.split('/')[-1]+"\t")

#fileLines = []
#for line in InputFile:
#	fileLines.append(line)

OutPutFile = open("JobInputDetails.dat","w")

for sample in Name.keys():
	#print "===>\t",sample
	InputFile = open("TotalEventInfo.dat","r")
	for line in InputFile:
		#print line.split()," -- ",sample
		if (line.split()[0]).find(sample) != -1:
			tmp2 = line
			#print "==> ",sample,"\t",line.split()
			for obj1 in OutputArray:
				if (obj1.split()[-1]).find(sample) != -1:
					tmp1 = obj1
					#if ('tmp1' in locals()) and ('tmp2' in locals()):
					#print tmp1.split(),"\t",tmp2.split()
					print '( '+tmp1.split()[-2],',\t"',tmp1.split()[-1],'",\t',tmp2.split()[-2],',\t',tmp2.split()[-1],')'
					OutPutFile.write('( '+str(tmp1.split()[-2])+',\t"'+str(tmp1.split()[-1])+'",\t'+str(tmp2.split()[-2])+',\t'+str(tmp2.split()[-1])+"),\n")
OutPutFile.close()
# cross-section	SampleDirName	nEvent	nNegEvent
"""
for sample in Name.keys():
	print "===>\t",sample
	for line in InputFile:
		if (line.split()[0]).find(sample) != -1:
			tmp2 = line
	for obj1 in OutputArray:
		if (obj1.split()[-1]).find(sample) != -1:
			tmp1 = obj1
	if ('tmp1' in locals()) and ('tmp2' in locals()):
		print tmp1.split(),"\t",tmp2.split()
		#print tmp1.split()[-2],"\t",tmp1.split()[-1],"\t",tmp2.split()[-2],'\t',tmp2.split()[-1]
# cross-section	SampleDirName	nEvent	nNegEvent
"""
