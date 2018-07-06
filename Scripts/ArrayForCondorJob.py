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

from pprint import pprint
#print "cross sections:" 
#pprint(Name)
#source = "/eos/uscms/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/"
source = "/eos/uscms/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples_New/"

#os.system('xrdfs root://cmseos.fnal.gov/ mkdir '+source+'HaddedFiles/')


ArrayDirPath = []
for root, dirs, filenames in os.walk(source):
	#print dirs
	for d in dirs:
		#print d
		ArrayDirPath.append(root+os.sep+d)

#print "=====	Print all dir	==========\n"
#pprint(ArrayDirPath)
print "Size of array = ",len(ArrayDirPath)

outScript = open("Summary_BaconNew.txt","w");

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
				nEvents = 1
				nNegEvents = 0
			else:
				print "==> ",directory,"\n\n"
				output = subprocess.check_output('python getnEvents.py -i  ' + directory + ' -t "Events" -o "test.txt"', shell=True) 
				print output
				nEvents = int(output.split()[3])
				nNegEvents = int(output.split()[6])
				#for root, dirs, filenames in os.walk(directory):
				#	for f in filenames:
				#		if f.endswith(".root"):
				#			tfile = ROOT.TFile(directory+os.sep+f)
				#			if not tfile or tfile.IsZombie():
				#				print "Failed or open tfile or IsZombie: ",directory+os.sep+f
				#				print ''
				#				exit(0)
				#			tree = tfile.Get('Events')
				#			if tree:
				#				if directory.find('Single') != -1:
				#					nEvents=1
				#				else:
				#					nEvents+=tree.GetEntries();
			print "==="
			print directory.split('/')[-1],"\t:\t  nEvents = ",nEvents,"\tnNegEvents = ",nNegEvents
			outScript.write(directory.split('/')[-1]+"\t:\t  nEvents = "+str(nEvents)+"\tnNegEvents = "+str(nNegEvents)+"\n\n")
			#print directory.split('/')[-1],"\t:\t  ",output.split()[1],output.split()[2],output.split()[3],"\t",output.split()[4],output.split()[5],output.split()[6]
			print "==="
			#outScript.write(directory.split('/')[-1]+"\t:\t  "+output.split()[1]+output.split()[2]+output.split()[3]+"\t"+output.split()[4]+output.split()[5]+output.split()[6]+"\n")
			nTotalEvents += nEvents
			nTotalNegEvents += nNegEvents
		#else:
		#	print "==> ",sample,"\t==\t Not found!!!"
	print "\t===> Total nEvents = ",nTotalEvents
	outScript.write("\t\t"+str(sample)+"\t\t Total nEvents = "+str(nTotalEvents)+"\t Tot neg events = "+str(nTotalNegEvents)+"\n\n")
	outScript.write('******\n\n')
	print "****"*20
outScript.close()
