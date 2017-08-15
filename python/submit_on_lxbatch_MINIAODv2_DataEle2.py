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

inputFolder = "/store/user/ksung/production/12a/";
outputFolder = currentDir+"/WWAnalysisRun2/output/";
exeName = currentDir+"/WWAnalysisRun2/produceWWNtuples"

dryRun = False;
doData = True;

category = ["el"];		# run this only for electrons



nameDataMu = [
"SingleMuonRun2016B_03Feb2017_ver1_v1",
"SingleMuonRun2016B_03Feb2017_ver2_v2",
"SingleMuonRun2016C_03Feb2017_v1",
"SingleMuonRun2016D_03Feb2017_v1",
"SingleMuonRun2016E_03Feb2017_v1",
"SingleMuonRun2016F_03Feb2017_v1",
"SingleMuonRun2016G_03Feb2017_v1",
"SingleMuonRun2016H_03Feb2017_ver2_v1",
"SingleMuonRun2016H_03Feb2017_ver3_v1"
];

nameDataEl = [
"SingleElectron_Run2016B-03Feb2017_ver1-v1",
"SingleElectron_Run2016B-03Feb2017_ver2-v2",
"SingleElectron_Run2016C-03Feb2017-v1",
"SingleElectron_Run2016D-03Feb2017-v1",
];


nameData = {"el": nameDataEl, "mu":nameDataMu};



for a in range(len(category)):
    
    if( doData ):
        for i in range(len(nameData[category[a]])):
            fn = "WWAnalysisRun2/Job/job_"+(nameData[category[a]])[i]+"_"+category[a];
            outScript = open(fn+".sh","w");
            command = "python "+currentDir+"/WWAnalysisRun2/python/produceWWNtuples.py -i "+inputFolder+" -n "+(nameData[category[a]])[i]+" -o WWTree_"+(nameData[category[a]])[i]+" -l "+category[a]+" -w 1. -no 1. --ismc 0 -trig 1";
            print command;
            outScript.write('#!/bin/bash');
            outScript.write("\n"+'cd '+CMSSWDir);
            outScript.write("\n"+'eval `scram runtime -sh`');
	    outScript.write("\n"+"cd WWAnalysis/WWAnalysisRun2");
            outScript.write("\n"+command);
	    outScript.write("\n"+"echo \"====> LISTING ALL FILES..... \"");
	    #outScript.write("\n"+"ls");
	    outScript.write("\n"+"echo \"====> MOVING OUTPUT ROOT FILES\"");
            outScript.write("\n"+"echo \"mv WWTree_"+(nameData[category[a]])[i]+".root "+outputFolder+"/output_"+category[a]+"\"");
            outScript.write("\n"+"mv WWTree_"+(nameData[category[a]])[i]+".root "+outputFolder+"/output_EleMu/");
	    outScript.write("\n"+"echo \"#########i########### Finishedi COPYING...\"");
            outScript.write("\n");
            outScript.close();
            os.system("chmod 777 "+currentDir+"/"+fn+".sh");
            command2 = "bsub -q 1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
            print command2
            if( dryRun != True ):
                os.system(command2);
                time.sleep(1)
