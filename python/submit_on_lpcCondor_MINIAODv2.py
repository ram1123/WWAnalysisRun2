#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess

currentDir = os.getcwd();
CMSSWDir =  currentDir+"/../";

#inputFolder = "/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/Feb142016";
inputFolder = "/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/";
outputFolder = currentDir+"/output/";
exePathName = currentDir+"/WWAnalysisRun2/produceWWNtuples"

lumi = 35900.0

dryRun = False;
doMC = True;
doData = False;

category = ["el"];
#category = ["el"];
#category = ["mu"];

samples = [
    #( 542.00,           "TTToSemilepton_powheg",                        0)
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo",	0)
    ]

nameDataMu = [
"data_mu_2016_runH_v2_9"
    ];

nameDataEl = [
"SingleElectron_Run2016B-03Feb2017_ver1-v1",
"SingleElectron_Run2016B-03Feb2017_ver2-v2",
"SingleElectron_Run2016C-03Feb2017-v1",
"SingleElectron_Run2016D-03Feb2017-v1"
];


inputlist="python/produceWWNtuples.py"

nameData = {"el": nameDataEl, "mu":nameDataMu};

#command = "python "+currentDir+"/WWAnalysisRun2/python/produceWWNtuples.py -i "+inputFolder+" $*";
command = "python produceWWNtuples.py -i "+inputFolder+" $*";

outScript = open("runstep2condor.sh","w");
outScript.write('#!/bin/bash');
outScript.write("\n"+'cd '+CMSSWDir);
outScript.write("\n"+'eval `scram runtime -sh`');
outScript.write("\n"+"cd -");
outScript.write("\n"+command);
outScript.write("\n");
outScript.close();
os.system("chmod 777 runstep2condor.sh");

outJDL = open("runstep2condor.jdl","w");
outJDL.write("Executable = runstep2condor.sh\n");
outJDL.write("Universe = vanilla\n");
outJDL.write("Requirements =FileSystemDomain==\"fnal.gov\" && Arch==\"X86_64\"");
outJDL.write("\n");
outJDL.write("Notification = ERROR\n");
outJDL.write("Should_Transfer_Files = YES\n");
outJDL.write("WhenToTransferOutput = ON_EXIT\n");
#outJDL.write("include : list-infiles.sh |\n");
outJDL.write("transfer_input_files = "+inputlist+"\n");
outJDL.write("x509userproxy = $ENV(X509_USER_PROXY)\n");

for a in range(len(category)):

    #MC
    if( doMC ):
        for i in range(len(samples)):
            outJDL.write("Output = "+str(samples[i][1])+"_"+category[a]+".stdout\n");
            outJDL.write("Error = "+str(samples[i][1])+"_"+category[a]+".stdout\n");
            outJDL.write("Arguments = -n "+str(samples[i][1])+" -o WWTree_"+str(samples[i][1])+"_"+category[a]+"_"+category[a]+" -w "+str(samples[i][0])+" -lumi "+str(lumi)+" --ismc 1 -trig 1 -c lpc\n");
            outJDL.write("Queue\n");
    
    #data
    if( doData ):
        for i in range(len(nameData[category[a]])):
            outJDL.write("Output = "+(nameData[category[a]])[i]+".stdout\n");
            outJDL.write("Error = "+(nameData[category[a]])[i]+".stdout\n");
            outJDL.write("Arguments = -n "+(nameData[category[a]])[i]+" -o WWTree_"+(nameData[category[a]])[i]+"_"+category[a]+" -w 1. -no 1. --ismc 0 -trig 1 -c lpc\n");
            outJDL.write("Queue\n");

outJDL.close();
print "===> Set Proxy Using:";
print "\tvoms-proxy-init --voms cms --valid 168:00";
print "\"condor_submit runstep2condor.jdl\" to submit";

