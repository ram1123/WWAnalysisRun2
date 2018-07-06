#!/usr/bin/python

import sys
import os
import subprocess


#source = "/eos/uscms/store/user/rasharma/MonteCarlo_Samples/MadGraph/EWK/"
#source = "/eos/uscms/store/user/rasharma/MonteCarlo_Samples/MadGraph/QCD/"
#source = "/eos/uscms/store/user/rasharma/MonteCarlo_Samples/MadGraph/EWKaQCD/"
source = "/eos/uscms/store/user/rasharma/MonteCarlo_Samples/MadGraph/EWK_Wlnujj_Ext/"

#os.system('xrdfs root://cmseos.fnal.gov/ mkdir '+source+'HaddedFiles/')


#SampleNames = [ "WPlepWMhadJJ_SM_LO_EWK_mjj100_pTj10", "WPhadWMlepJJ_SM_LO_EWK_mjj100_pTj10",
#		"WPlepWPhadJJ_SM_LO_EWK_mjj100_pTj10", "WMlepWMhadJJ_SM_LO_EWK_mjj100_pTj10",
#		"WPlepZhadJJ_SM_LO_EWK_mjj100_pTj10",  "WMlepZhadJJ_SM_LO_EWK_mjj100_pTj10",
#		"WPhadZlepJJ_SM_LO_EWK_mjj100_pTj10", "WMhadZlepJJ_SM_LO_EWK_mjj100_pTj10",
#		"ZlepZhadJJ_SM_LO_EWK_mjj100_pTj10"
#		]
#SampleNames = [ "WPlepWMhadJJ_SM_LO_QCD_mjj100_pTj10", "WPhadWMlepJJ_SM_LO_QCD_mjj100_pTj10",
#		"WPlepWPhadJJ_SM_LO_QCD_mjj100_pTj10", "WMlepWMhadJJ_SM_LO_QCD_mjj100_pTj10",
#		"WPlepZhadJJ_SM_LO_QCD_mjj100_pTj10",  "WMlepZhadJJ_SM_LO_QCD_mjj100_pTj10",
#		"WPhadZlepJJ_SM_LO_QCD_mjj100_pTj10", "WMhadZlepJJ_SM_LO_QCD_mjj100_pTj10",
#		"ZlepZhadJJ_SM_LO_QCD_mjj100_pTj10"
#		]
#SampleNames = [ "WPlepWMhadJJ_SM_LO_EWKaQCD_mjj100_pTj10", "WPhadWMlepJJ_SM_LO_EWKaQCD_mjj100_pTj10",
#		"WPlepWPhadJJ_SM_LO_EWKaQCD_mjj100_pTj10", "WMlepWMhadJJ_SM_LO_EWKaQCD_mjj100_pTj10",
#		"WPlepZhadJJ_SM_LO_EWKaQCD_mjj100_pTj10",  "WMlepZhadJJ_SM_LO_EWKaQCD_mjj100_pTj10",
#		"WPhadZlepJJ_SM_LO_EWKaQCD_mjj100_pTj10", "WMhadZlepJJ_SM_LO_EWKaQCD_mjj100_pTj10",
#		"ZlepZhadJJ_SM_LO_EWKaQCD_mjj100_pTj10"
#		]
SampleNames = [ 
		"WMhadElenuJJ_SM_LO_EWK_mjj100_pTj10",  "WMhadTaunuJJ_SM_LO_EWK_mjj100_pTj10",
		"WPhadMunuJJ_SM_LO_EWK_mjj100_pTj10",	"WMhadMunuJJ_SM_LO_EWK_mjj100_pTj10",
		"WPhadElenuJJ_SM_LO_EWK_mjj100_pTj10",	"WPhadTaunuJJ_SM_LO_EWK_mjj100_pTj10"
		]

ArrayDirPath = []
Count = 0
for root, dirs, filenames in os.walk(source):
	for files in filenames:
		for samples in SampleNames:
			if samples in root:
				#if files.endswith(".txt"):
				#	command = 'xrdcp '+root + os.sep + files + " EWK_Wlnujj_Samples/"+samples + os.sep + os.path.splitext(files)[0]+"_"+str(Count)+os.path.splitext(files)[1]
				#	print command
				#	os.system(command)
				#	print "\n\n"
				#	Count += 1
				if files.endswith(".gz"):
					command = 'xrdcp '+root + os.sep + files + " EWK_Wlnujj_Samples/"+samples + os.sep + files.split(".")[0] + "_" + str(Count) +".lhe.gz" 
					print command
					os.system(command)
					print "\n\n"
					Count += 1


### total 526 .txt files
