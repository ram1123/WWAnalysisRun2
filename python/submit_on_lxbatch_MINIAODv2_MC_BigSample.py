#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess

currentDir = os.getcwd();
CMSSWDir = currentDir+"/../";

lumi = 35900.0

inputFolder = "/store/cmst3/group/monojet/production/12a/";	# Path of data files
outputFolder = currentDir+"/WWAnalysisRun2/output/";
exeName = currentDir+"/WWAnalysisRun2/produceWWNtuples"

dryRun = False;
doMC = False;
doData = True;	# For data run another script

#category = ["EleMu"];
#category = ["el"];
category = ["mu"];



inputFolder = "/store/cmst3/group/monojet/production/12/";	# Path of MC files
samples = [
    ###( 27990000,	"QCD_HT100to200_13TeV",		0),
    ###( 1712000,	"QCD_HT200to300_13TeV",		0),
    #( 1712000,	"QCD_HT200to300_13TeV_ext",	0),
    ( 347700,	"QCD_HT300to500_13TeV",		0),
    #( 347700,	"QCD_HT300to500_13TeV_ext",	0),
    #( 32100,	"QCD_HT500to700_13TeV",		0),
    #( 32100,	"QCD_HT500to700_13TeV_ext",	0),
    ###( 6831,	"QCD_HT700to1000_13TeV",	0),
    #( 6831,	"QCD_HT700to1000_13TeV_ext",	0),
    #( 1207,	"QCD_HT1000to1500_13TeV",	0),
    #( 1207,	"QCD_HT1000to1500_13TeV_ext",	0),
    #( 119.9,	"QCD_HT1500to2000_13TeV",	0),
    #( 119.9,	"QCD_HT1500to2000_13TeV_ext",	0),
    #( 25.24,	"QCD_HT2000toInf_13TeV",	0),
    #( 25.24,	"QCD_HT2000toInf_13TeV_ext",	0),
    #( 49.997,		"WWToLNuQQ_13TeV_powheg",			0),
    ###( 49.997,		"WWToLNuQQ_13TeV_powheg_ext",			0),
    #( 61526.7,		"WJetsToLNu_13TeV",				0),
    ##( 			"WJetsToLNu_HT_70To100_13TeV"
    #( 1506.4,		"WJetsToLNu_HT_100To200_13TeV", 		0),
    #( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext1",		0),
    ###( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext2",		0),
    #( 435.237,		"WJetsToLNu_HT_200To400_13TeV",			0),
    #( 435.237,		"WJetsToLNu_HT_200To400_13TeV_ext1",		0),
    ###( 435.237,		"WJetsToLNu_HT_200To400_13TeV_ext2",		0),
    #( 59.1811,		"WJetsToLNu_HT_400To600_13TeV",			0),
    ###( 59.1811,		"WJetsToLNu_HT_400To600_13TeV_ext1",		0),
    #( 14.5805,		"WJetsToLNu_HT_600To800_13TeV",			0),
    #( 14.5805,		"WJetsToLNu_HT_600To800_13TeV_ext1",		0),
    #( 6.65621,		"WJetsToLNu_HT_800To1200_13TeV",		0),
    ##( 6.65621,		"WJetsToLNu_HT_800To1200_13TeV_ext1",	
    #( 1.60809,		"WJetsToLNu_HT_1200To2500_13TeV",		0),
    ###( 1.60809,		"WJetsToLNu_HT_1200To2500_13TeV_ext1",		0),
    #( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV",		0),
    #( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV_ext1",		0),
    #( 49.997,		"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",0),
    #( 49.997,		"WWToLNuQQ_13TeV_powheg",			0),
    #( 49.997,		"WWToLNuQQ_13TeV_powheg_ext",			0),
    ###( 10.71,		"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",0),
    ###( 3.22,		"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",	0),
    #( 542.0,		"TTToSemilepton_powheg",			0),
    #( 11.36,		"ST_s_channel_4f_leptonDecays_13TeV_amcatnlo_pythia8_TuneCUETP8M1",	0),
    #( 80.95,		"ST_t_channel_antitop_4f_inclusiveDecays_13TeV_powhegV2_madspin_pythia8_TuneCUETP8M1",	0),
    ##( 136.02,		"ST_t_channel_top_4f_inclusiveDecays_13TeV_powhegV2_madspin_pythia8_TuneCUETP8M1",	0),
    #( 19.46,		"ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",	0),
    #( 19.46,		"ST_tW_top_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",		0),
    #( 136.02,		"ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV_powhegV2_madspin",	0)
    ]

nameDataMu = [
#"SingleMuonRun2016B_03Feb2017_ver1_v1",
#"SingleMuonRun2016B_03Feb2017_ver2_v2",
#"SingleMuonRun2016C_03Feb2017_v1",
#"SingleMuonRun2016D_03Feb2017_v1",
#"SingleMuonRun2016E_03Feb2017_v1",
"SingleMuonRun2016F_03Feb2017_v1"
#"SingleMuonRun2016G_03Feb2017_v1",
#"SingleMuonRun2016H_03Feb2017_ver2_v1",
#"SingleMuonRun2016H_03Feb2017_ver3_v1"
];

nameDataEl = [
#"SingleElectronRun2016B_23Sep2016_v2",
#"SingleElectronRun2016C_23Sep2016_v1",
#"SingleElectronRun2016H_PromptReco_v1",
"SingleElectronRun2016H_PromptReco_v3_test"
];


nameData = {"el": nameDataEl, "mu":nameDataMu};

#import commands
import subprocess

#MC
if( doMC ):
    for i in range(len(samples)):
        print samples[i][1]
        command = "eos find -f "+inputFolder+"/"+samples[i][1]
        print command
        count=0
        with os.popen(command) as pipe:
        	for line in pipe:
    		    if line.strip().find("log")==-1:
    			if line.strip().find("fail")==-1:
    				count+=1
				#print count,"==> ",line.strip()
        			fn = "WWAnalysisRun2/Job_new/job_"+samples[i][1]+"_"+str(count);
        			outScript = open(fn+".sh","w");
        			command = "python "+currentDir+"/WWAnalysisRun2/python/produceWWNtuples.py -i "+str(line.strip())+" -o WWTree_"+str(samples[i][1])+"_"+str(count)+" -w "+str(samples[i][0])+" -lumi "+str(lumi)+" --ismc 1 -trig 0 -n tr ";
        			#command = "python "+currentDir+"/WWAnalysisRun2/python/produceWWNtuples.py -i "+str(line.strip())+" -o WWTree_"+str(samples[i][1])+"_"+category[a]+" -l "+category[a]+" -w "+str(samples[i][0])+" -lumi "+str(lumi)+" --ismc 1 -trig 0";
        			print command;
        			outScript.write('#!/bin/bash');
        			outScript.write("\n"+'cd '+CMSSWDir);
        			outScript.write("\n"+'eval `scram runtime -sh`');
        			outScript.write("\n"+"cd WWAnalysis/WWAnalysisRun2");
        			outScript.write("\n"+command);
        			outScript.write("\n"+"echo \"====> LISTING ALL FILES..... \"");
        			outScript.write("\n"+"ls -ltrh");
        			outScript.write("\n"+"echo \"====> MOVING OUTPUT ROOT FILES\"");
				#outScript.write("\n"+"mkdir "+samples[i][1]);
        			outScript.write("\n"+"echo \"mv WWTree_"+(samples[i][1])+"_"+str(count)+".root "+outputFolder+"/output_"+category[0]+"/"+samples[i][1]+"\"");
        			outScript.write("\n"+"mv WWTree_"+(samples[i][1])+"_"+str(count)+".root "+outputFolder+"/output_"+category[0]+"/"+samples[i][1]);
        			outScript.write("\n"+"echo \"#########i########### Finished COPYING...\"");
        			outScript.write("\n");
        			os.system("chmod 777 "+currentDir+"/"+fn+".sh");
        			command2 = "bsub -q 1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh -J "+str(samples[i][1][0:9]);
				outScript.write("\n"+"Bjob submission command: "+command2);
        			print command2
        			outScript.close();
        			if( dryRun != True ):
				    print "\n======= JOB: ",count," ======="
				    dirCreate="mkdir "+outputFolder+"/output_"+category[0]+"/"+samples[i][1]
				    if count==1:
				    	os.system(dirCreate);
				    os.system(command2);
        			    #time.sleep(3)

#data
if( doData ):
    for a in range(len(category)):
        for i in range(len(nameData[category[a]])):
            command = "eos find -f "+inputFolder+"/"+(nameData[category[a]])[i]
	    print command
	    count=0
	    with os.popen(command) as pipe:
	          for line in pipe:
		     if line.strip().find("log")==-1:
		        if line.strip().find("fail")==-1:
			   count+=1
			   fn = "WWAnalysisRun2/Job_new/job_"+(nameData[category[a]])[i]+"_"+str(count);
        		   outScript = open(fn+".sh","w");
            		   command = "python "+currentDir+"/WWAnalysisRun2/python/produceWWNtuples.py -i "+str(line.strip())+" -o WWTree_"+(nameData[category[a]])[i]+"_"+str(count)+" -l "+category[a]+" -w 1. -no 1. -lumi "+str(lumi)+" --ismc 0 -trig 1 -n tr";
            		   print command;
            		   outScript.write('#!/bin/bash');
            		   outScript.write("\n"+'cd '+CMSSWDir);
            		   outScript.write("\n"+'eval `scram runtime -sh`');
	    		   outScript.write("\n"+"cd WWAnalysis/WWAnalysisRun2");
            		   outScript.write("\n"+command);
	    		   outScript.write("\n"+"echo \"====> LISTING ALL FILES..... \"");
	    		   #outScript.write("\n"+"ls");
	    		   outScript.write("\n"+"echo \"====> MOVING OUTPUT ROOT FILES\"");
            		   outScript.write("\n"+"echo \"mv WWTree_"+(nameData[category[a]])[i]+"_"+str(count)+".root "+outputFolder+"/"+(nameData[category[a]])[i]+"\"");
            		   outScript.write("\n"+"mv WWTree_"+(nameData[category[a]])[i]+"_"+str(count)+".root "+outputFolder+"/"+(nameData[category[a]])[i]);
	    		   outScript.write("\n"+"echo \"#########i########### Finishedi COPYING...\"");
            		   outScript.write("\n");
            		   outScript.close();
            		   os.system("chmod 777 "+currentDir+"/"+fn+".sh");
            		   command2 = "bsub -q 1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
            		   print command2
            		   if( dryRun != True ):
			       print "\n======= JOB: ",count," ======="
			       if count==1:
			          dirCreate="mkdir "+outputFolder+"/"+(nameData[category[a]])[i]
			       	  os.system(dirCreate);
			       os.system(command2);
