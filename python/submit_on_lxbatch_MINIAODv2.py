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

inputFolder = "/store/cmst3/group/monojet/production/12/";	# Path of MC files
#inputFolder = "/store/cmst3/group/monojet/production/12a/";	# Path of data files
outputFolder = currentDir+"/WWAnalysisRun2/output/";
exeName = currentDir+"/WWAnalysisRun2/produceWWNtuples"

dryRun = False;
doMC = False;
doData = True;

category = ["mu","el"];
#category = ["el"];
#category = ["mu"];



samples = [
    ( 49.997,		"WWToLNuQQ_13TeV_powheg",			0),
    ( 49.997,		"WWToLNuQQ_13TeV_powheg_ext",			0),
    ( 61526.7,		"WJetsToLNu_13TeV",				0),
    #( 			"WJetsToLNu_HT_70To100_13TeV"
    ( 1506.4,		"WJetsToLNu_HT_100To200_13TeV", 		0),
    ( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext1",		0),
    ( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext2",		0),
    ( 435.237,		"WJetsToLNu_HT_200To400_13TeV",			0),
    ( 435.237,		"WJetsToLNu_HT_200To400_13TeV_ext1",		0),
    ( 435.237,		"WJetsToLNu_HT_200To400_13TeV_ext2",		0),
    ( 59.1811,		"WJetsToLNu_HT_400To600_13TeV",			0),
    ( 59.1811,		"WJetsToLNu_HT_400To600_13TeV_ext1",		0),
    ( 14.5805,		"WJetsToLNu_HT_600To800_13TeV",			0),
    ( 14.5805,		"WJetsToLNu_HT_600To800_13TeV_ext1",		0),
    ( 6.65621,		"WJetsToLNu_HT_800To1200_13TeV",		0),
    #( 6.65621,		"WJetsToLNu_HT_800To1200_13TeV_ext1",	
    ( 1.60809,		"WJetsToLNu_HT_1200To2500_13TeV",		0),
    ( 1.60809,		"WJetsToLNu_HT_1200To2500_13TeV_ext1",		0),
    ( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV",		0),
    ( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV_ext1",		0),
    ( 49.997,		"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",0),
    ( 49.997,		"WWToLNuQQ_13TeV_powheg",			0),
    ( 49.997,		"WWToLNuQQ_13TeV_powheg_ext",			0),
    ( 10.71,		"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",0),
    ( 3.22,		"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",	0),
    ( 831.76,		"TTToSemilepton_powheg",			0),
    ( 11.36,		"ST_s_channel_4f_leptonDecays_13TeV_amcatnlo_pythia8_TuneCUETP8M1",	0),
    ( 80.95,		"ST_t_channel_antitop_4f_inclusiveDecays_13TeV_powhegV2_madspin_pythia8_TuneCUETP8M1",	0),
    ( 136.02,		"ST_t_channel_top_4f_inclusiveDecays_13TeV_powhegV2_madspin_pythia8_TuneCUETP8M1",	0),
    ( 19.46,		"ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",	0),
    ( 19.46,		"ST_tW_top_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",		0)
    ]

nameDataMu = [
"SingleMuonRun2016B_23Sep2016_v1",
"SingleMuonRun2016B_23Sep2016_v3",
"SingleMuonRun2016C_23Sep2016_v1",
"SingleMuonRun2016D_23Sep2016_v1",
"SingleMuonRun2016E_23Sep2016_v1",
"SingleMuonRun2016F_23Sep2016_v1",
"SingleMuonRun2016G_23Sep2016_v1",
"SingleMuonRun2016H_PromptReco_v1",
"SingleMuonRun2016H_PromptReco_v2",
"SingleMuonRun2016H_PromptReco_v3"

#"SingleMuonRun2016B_03Feb2017_ver1_v1",
#"SingleMuonRun2016B_03Feb2017_ver2_v2",
#"SingleMuonRun2016C_03Feb2017_v1",
#"SingleMuonRun2016D_03Feb2017_v1",
#"SingleMuonRun2016E_03Feb2017_v1",
#"SingleMuonRun2016F_03Feb2017_v1",
#"SingleMuonRun2016G_03Feb2017_v1",
#"SingleMuonRun2016H_03Feb2017_ver2_v1",
#"SingleMuonRun2016H_03Feb2017_ver3_v1"
];

nameDataEl = [
"SingleElectronRun2016B_23Sep2016_v2",
"SingleElectronRun2016C_23Sep2016_v1",
"SingleElectronRun2016H_PromptReco_v1",
"SingleElectronRun2016H_PromptReco_v3",

#"SingleElectron_Run2016C-03Feb2017-v1",
#"SingleElectron_Run2016E-03Feb2017-v1",
#"SingleElectronRun2016F_03Feb2017_v1",
#"SingleElectron_Run2016G-03Feb2017-v1",
#"SingleElectronRun2016H_03Feb2017_ver2_v1",
#"SingleElectronRun2016H_03Feb2017_ver3_v1",
#"SingleElectronRun2016F_03Feb2017_v1",
#"SingleElectronRun2016H_03Feb2017_ver2_v1",
#"SingleElectronRun2016H_03Feb2017_ver3_v1"
];


nameData = {"el": nameDataEl, "mu":nameDataMu};



for a in range(len(category)):
    
    #MC
    if( doMC ):
        for i in range(len(samples)):
            fn = "WWAnalysisRun2/Job/job_"+samples[i][1]+"_"+category[a];
            outScript = open(fn+".sh","w");
            command = "python "+currentDir+"/WWAnalysisRun2/python/produceWWNtuples.py -i "+inputFolder+" -n "+str(samples[i][1])+" -o WWTree_"+str(samples[i][1])+" -l "+category[a]+" -w "+str(samples[i][0])+" -mass "+str(samples[i][2])+" --ismc 1 -trig 0";
            print command;
            outScript.write('#!/bin/bash');
            outScript.write("\n"+'cd '+CMSSWDir);
            outScript.write("\n"+'eval `scram runtime -sh`');
            outScript.write("\n"+"cd WWAnalysis/WWAnalysisRun2");
            outScript.write("\n"+command);
	    outScript.write("\n"+"echo \"====> LISTING ALL FILES..... \"");
	    outScript.write("\n"+"ls");
	    outScript.write("\n"+"echo \"====> MOVING OUTPUT ROOT FILES\"");
	    outScript.write("\n"+"cp WWTree_"+(samples[i][1])+"_"+category[a]+".root "+outputFolder+"/output_"+category[a]);
	    outScript.write("\n"+"echo \"#########i########### Finished COPYING...\"");
            outScript.write("\n");
            outScript.close();
            os.system("chmod 777 "+currentDir+"/"+fn+".sh");
            command2 = "bsub -q 1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
            print command2
            if( dryRun != True ):
                os.system(command2);
                time.sleep(3)
    
    #data
    if( doData ):
        for i in range(len(nameData[category[a]])):
            fn = "WWAnalysisRun2/Job/job_"+(nameData[category[a]])[i]+"_"+category[a];
            outScript = open(fn+".sh","w");
            command = "python "+currentDir+"/WWAnalysisRun2/python/produceWWNtuples.py -i "+inputFolder+" -n "+(nameData[category[a]])[i]+" -o WWTree_"+(nameData[category[a]])[i]+" -l "+category[a]+" -w 1. -no 1. -mass 0 --ismc 0 -trig 1";
            print command;
            outScript.write('#!/bin/bash');
            outScript.write("\n"+'cd '+CMSSWDir);
            outScript.write("\n"+'eval `scram runtime -sh`');
	    outScript.write("\n"+"cd WWAnalysis/WWAnalysisRun2");
            outScript.write("\n"+command);
	    outScript.write("\n"+"echo \"====> LISTING ALL FILES..... \"");
	    outScript.write("\n"+"ls");
	    outScript.write("\n"+"echo \"====> MOVING OUTPUT ROOT FILES\"");
            outScript.write("\n"+"cp WWTree_"+(nameData[category[a]])[i]+".root "+outputFolder+"/output_"+category[a]);
	    outScript.write("\n"+"echo \"#########i########### Finishedi COPYING...\"");
            outScript.write("\n");
            outScript.close();
            os.system("chmod 777 "+currentDir+"/"+fn+".sh");
            command2 = "bsub -q 1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
            print command2
            if( dryRun != True ):
                os.system(command2);
                time.sleep(1)
