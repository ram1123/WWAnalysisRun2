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

#inputFolder = "/store/cmst3/group/monojet/production/12";
#inputFolder = "/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples_New/";
#inputFolder = "/eos/cms/store/user/ksung/production/12/"
inputFolder = "/eos/uscms/store/user/aapyan2/Run2/"
outputFolder = "/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/";

dryRun = False;
doMC = False;
doData = True;

sampleName = [
"ChargedHiggsToWZToLLQQ_M300_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLLQQ_M200_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLLQQ_M400_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLLQQ_M600_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLLQQ_M800_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLNuQQ_M300_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLNuQQ_M500_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLNuQQ_M700_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLNuQQ_M200_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLNuQQ_M800_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLNuQQ_M400_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLNuQQ_M2000_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLLQQ_M700_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLLQQ_M500_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLLQQ_M1500_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLNuQQ_M1500_13TeV-madgraph-pythia8",
"ChargedHiggsToWZToLNuQQ_M600_13TeV-madgraph-pythia8"
#'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1',
#'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2',
#'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_3',
#'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_4',
#'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_5',
#'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_6',
#'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_7',
#'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_8',
#'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_9',
#'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1',
#'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2',
#'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_3',
#'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_4',
#'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_5',
#'DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1',
#'DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2',
#'DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
#'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
#'DYToLL_0J_13TeV-amcatnloFXFX-pythia8',
#'DYToLL_1J_13TeV-amcatnloFXFX-pythia8_1',
#'DYToLL_1J_13TeV-amcatnloFXFX-pythia8_2',
#'DYToLL_1J_13TeV-amcatnloFXFX-pythia8_3',
#'DYToLL_1J_13TeV-amcatnloFXFX-pythia8_4',
#'DYToLL_1J_13TeV-amcatnloFXFX-pythia8_5',
#'DYToLL_1J_13TeV-amcatnloFXFX-pythia8_6',
#'DYToLL_1J_13TeV-amcatnloFXFX-pythia8_7',
#'DYToLL_2J_13TeV-amcatnloFXFX-pythia8_1',
#'DYToLL_2J_13TeV-amcatnloFXFX-pythia8_2',
#'DYToLL_2J_13TeV-amcatnloFXFX-pythia8_3',
#'DYToLL_2J_13TeV-amcatnloFXFX-pythia8_4',
#'DYToLL_2J_13TeV-amcatnloFXFX-pythia8_5',
#'WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WminusTo2JZTo2LJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WminusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WminusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WminusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WminusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WplusTo2JWminusToLNuJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WplusTo2JWminusToLNuJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WplusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WplusTo2JZTo2LJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WplusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WplusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WplusToLNuWplusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WplusToLNuWplusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WplusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'WplusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'ZTo2LZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8',
#'ZTo2LZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8'
]

command1 = "xrdfs root://cmseos.fnal.gov/ mkdir "+outputFolder+"/$1/"
command2 = "xrdcp -r root://cmseos.fnal.gov/"+inputFolder+"/$1/  root://cmseos.fnal.gov/"+outputFolder+"/$1/"
#command2 = "xrdcp -r -f -s root://eoscms.cern.ch/"+inputFolder+"/$1/  root://cmseos.fnal.gov/"+outputFolder+"/$1/"

outScript = open("runCopycondor.sh","w");
outScript.write('#!/bin/bash');
outScript.write("\n"+'cd '+CMSSWDir);
outScript.write("\n"+'eval `scram runtime -sh`');
outScript.write("\n"+"cd -");
outScript.write("\n"+"echo \"sample Name : \"$1");
outScript.write("\n"+command1);
outScript.write("\n"+command2);
outScript.write("\n");
outScript.close();
os.system("chmod 777 runCopycondor.sh");

outJDL = open("runCopycondor.jdl","w");
outJDL.write("Executable = runCopycondor.sh\n");
outJDL.write("Universe = vanilla\n");
outJDL.write("Requirements =FileSystemDomain==\"fnal.gov\" && Arch==\"X86_64\"");
outJDL.write("\n");
outJDL.write("Notification = ERROR\n");
outJDL.write("Should_Transfer_Files = YES\n");
outJDL.write("WhenToTransferOutput = ON_EXIT\n");
#outJDL.write("include : list-infiles.sh |\n");
outJDL.write("x509userproxy = $ENV(X509_USER_PROXY)\n");

for i in range(len(sampleName)):
    outJDL.write("Output = "+str(sampleName[i])+"_Copy.stdout\n");
    outJDL.write("Error = "+str(sampleName[i])+"_Copy.stdout\n");
    outJDL.write("Arguments = "+str(sampleName[i])+"\n");
    outJDL.write("Queue\n");

outJDL.close();
print "===> Set Proxy Using:";
print "\tvoms-proxy-init --voms cms --valid 168:00";
print "\"condor_submit runCopycondor.jdl\" to submit";

