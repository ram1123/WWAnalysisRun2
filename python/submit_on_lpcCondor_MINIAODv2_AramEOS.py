#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess
import tarfile
import datetime
import commands

currentDir = os.getcwd();
CMSSWDir =  currentDir+"/../";

inputFolder = "/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples_New/"
TestRun = 0

doMC = True;
doData = False;
category = ["el","mu"];

lumi = 35900.0

# Get date and time for output directory
## ADD "test" IN OUTPUT FOLDER IF YOU ARE TESTING SO THAT LATER YOU REMEMBER TO WHICH DIRECTORY YOU HAVE TO REMOVE FROM EOS
if TestRun:
	outputFolder = "/store/user/rasharma/SecondStep/WWTree_"+datetime.datetime.now().strftime('%Y-%m-%d_%Hh%M')+"_TEST/";
	OutputLogPath = "OutPut_Logs/Logs_" + datetime.datetime.now().strftime('%Y-%m-%d_%Hh%M') + "_TEST";
else:
	outputFolder = "/store/user/rasharma/SecondStep/WWTree_"+datetime.datetime.now().strftime('%Y-%m-%d_%Hh%M');
	OutputLogPath = "OutPut_Logs/Logs_" + datetime.datetime.now().strftime('%Y-%m-%d_%Hh%M');

outputFolder = "/store/user/rasharma/SecondStep/WWTree_2017_12_06_11h55/" 
OutputLogPath = "OutPut_Logs/Logs_2017_12_06_11h55/"

print "Name of output dir: ",outputFolder
# create a directory on eos
#os.system('xrdfs root://cmseos.fnal.gov/ mkdir ' + outputFolder)
# create directory in pwd for log files
#os.system('mkdir -p ' + OutputLogPath)

# Function to create a tar file
def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))

# Get CMSSW directory path and name
cmsswDirPath = commands.getstatusoutput('echo ${CMSSW_BASE}')
CMSSWRel = os.path.basename(cmsswDirPath[1])

print "CMSSW release used : ",CMSSWRel

# create tarball of present working CMSSW base directory
os.system('rm CMSSW*.tgz')
make_tarfile(CMSSWRel+".tgz", cmsswDirPath[1])

# send the created tarball to eos
os.system('xrdcp -f ' + CMSSWRel+".tgz" + ' root://cmseos.fnal.gov/' + outputFolder + '/' + CMSSWRel+".tgz")
os.system('xrdcp -f ThingsUpdated.txt root://cmseos.fnal.gov/' + outputFolder)
os.system('cp ThingsUpdated.txt ' + OutputLogPath)

samples = [
    ( 0.9114,	"WplusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	1991227,   0),
    ( 0.9107,	"WplusTo2JWminusToLNuJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	1983847,   0),
    ( 0.0879,	"WplusToLNuWplusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	198848,   0),
    ( 0.0326,	"WminusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	189560,   0),
    ( 0.1825,	"WplusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		393171,   0),
    ( 0.0540,	"WplusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		198922,   0),
    ( 0.1000,	"WminusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		199542,   0),
    ( 0.0298,	"WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		189086,   0),
    ( 0.0159,	"ZTo2LZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		99997,   0),
    ( 5.546,	"WplusTo2JWminusToLNuJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	3994663,   0),
    ( 5.558,	"WplusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	3949170,   0),
    ( 0.086,	"WplusToLNuWplusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	99992,   0),
    ( 0.038,	"WminusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	99657,   0),
    ( 2.159,	"WplusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		1991348,   0),
    ( 0.640,	"WplusTo2JZTo2LJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		499432,   0),
    ( 1.302,	"WminusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		981540,   0),
    ( 0.387,	"WminusTo2JZTo2LJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		489280,   0),
    ( 0.376,	"ZTo2LZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		49999,   0),
    ( 4274.1645,	"DYToLL_0J_13TeV-amcatnloFXFX-pythia8",			1749590,	161090),
    ( 1115.0463,	"DYToLL_1J_13TeV-amcatnloFXFX-pythia8_1",		35950579,	9808428),
    ( 1115.0463,	"DYToLL_1J_13TeV-amcatnloFXFX-pythia8_2",		35950579,	9808428),
    ( 1115.0463,	"DYToLL_1J_13TeV-amcatnloFXFX-pythia8_3",		35950579,	9808428),
    ( 1115.0463,	"DYToLL_1J_13TeV-amcatnloFXFX-pythia8_4",		35950579,	9808428),
    ( 1115.0463,	"DYToLL_1J_13TeV-amcatnloFXFX-pythia8_5",		35950579,	9808428),
    ( 1115.0463,	"DYToLL_1J_13TeV-amcatnloFXFX-pythia8_6",		35950579,	9808428),
    ( 1115.0463,	"DYToLL_1J_13TeV-amcatnloFXFX-pythia8_7",		35950579,	9808428),
    ( 222.57608,	"DYToLL_2J_13TeV-amcatnloFXFX-pythia8_1",		21571879,	7649488),
    ( 222.57608,	"DYToLL_2J_13TeV-amcatnloFXFX-pythia8_2",		21571879,	7649488),
    ( 222.57608,	"DYToLL_2J_13TeV-amcatnloFXFX-pythia8_3",		21571879,	7649488),
    ( 222.57608,	"DYToLL_2J_13TeV-amcatnloFXFX-pythia8_4",		21571879,	7649488),
    ( 222.57608,	"DYToLL_2J_13TeV-amcatnloFXFX-pythia8_5",		21571879,	7649488),
    ( 45.446723,	"DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",	2798791,	0),
    ( 97.554176,	"DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1",	4866978,	0),
    ( 97.554176,	"DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2",	4866978,	0),
    ( 314.57229,	"DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1",	19296117,	0),
    ( 314.57229,	"DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2",	19296117,	0),
    ( 314.57229,	"DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_3",	19296117,	0),
    ( 314.57229,	"DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_4",	19296117,	0),
    ( 314.57229,	"DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_5",	19296117,	0),
    ( 1103.0878,	"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1",	52602172,	0),
    ( 1103.0878,	"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2",	52602172,	0),
    ( 1103.0878,	"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_3",	52602172,	0),
    ( 1103.0878,	"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_4",	52602172,	0),
    ( 1103.0878,	"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_5",	52602172,	0),
    ( 1103.0878,	"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_6",	52602172,	0),
    ( 1103.0878,	"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_7",	52602172,	0),
    ( 1103.0878,	"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_8",	52602172,	0),
    ( 1103.0878,	"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_9",	52602172,	0),
    #( 4435.5258,	"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",	0)
        ]

nameDataMu = [
    "SingleMuonRun2016B_03Feb2017_ver1_v1",
    "SingleMuonRun2016B_03Feb2017_ver2_v2_1",
    "SingleMuonRun2016B_03Feb2017_ver2_v2_2",
    "SingleMuonRun2016B_03Feb2017_ver2_v2_3",
    "SingleMuonRun2016C_03Feb2017_v1",
    "SingleMuonRun2016D_03Feb2017_v1_1",
    "SingleMuonRun2016D_03Feb2017_v1_2",
    "SingleMuonRun2016E_03Feb2017_v1_1",
    "SingleMuonRun2016E_03Feb2017_v1_2",
    "SingleMuonRun2016F_03Feb2017_v1_1",
    "SingleMuonRun2016F_03Feb2017_v1_2",
    "SingleMuonRun2016G_03Feb2017_v1_1",
    "SingleMuonRun2016G_03Feb2017_v1_2",
    "SingleMuonRun2016G_03Feb2017_v1_3",
    "SingleMuonRun2016H_03Feb2017_ver2_v1_1",
    "SingleMuonRun2016H_03Feb2017_ver2_v1_2",
    "SingleMuonRun2016H_03Feb2017_ver2_v1_3",
    "SingleMuonRun2016H_03Feb2017_ver3_v1"
    ];

nameDataEl = [
"SingleElectron_Run2016B-03Feb2017_ver1-v1",
"SingleElectron_Run2016B-03Feb2017_ver2-v2_1",
"SingleElectron_Run2016B-03Feb2017_ver2-v2_2",
"SingleElectron_Run2016B-03Feb2017_ver2-v2_3",
"SingleElectron_Run2016B-03Feb2017_ver2-v2_4",
"SingleElectron_Run2016B-03Feb2017_ver2-v2_5",
"SingleElectron_Run2016C-03Feb2017-v1",
"SingleElectron_Run2016D-03Feb2017-v1_1",
"SingleElectron_Run2016D-03Feb2017-v1_2",
"SingleElectron_Run2016D-03Feb2017-v1_3",
"SingleElectron_Run2016E-03Feb2017-v1_1",
"SingleElectron_Run2016E-03Feb2017-v1_2",
"SingleElectronRun2016F_03Feb2017_v1",
"SingleElectron_Run2016G-03Feb2017-v1_1",
"SingleElectron_Run2016G-03Feb2017-v1_2",
"SingleElectron_Run2016G-03Feb2017-v1_3",
"SingleElectronRun2016H_03Feb2017_ver2_v1_1",
"SingleElectronRun2016H_03Feb2017_ver2_v1_2",
"SingleElectronRun2016H_03Feb2017_ver3_v1"
];


inputlist = "runstep2condor.sh, python/produceWWNtuples.py"

nameData = {"el": nameDataEl, "mu":nameDataMu};

command = "python python/produceWWNtuples.py -i "+inputFolder+" $*";

outScript = open("runstep2condor.sh","w");
outScript.write('#!/bin/bash');
outScript.write("\n"+'echo "Starting job on " `date`');
outScript.write("\n"+'echo "Running on: `uname -a`"');
outScript.write("\n"+'echo "System software: `cat /etc/redhat-release`"');
outScript.write("\n"+'source /cvmfs/cms.cern.ch/cmsset_default.sh');
outScript.write("\n"+'### copy the input root files if they are needed only if you require local reading');
outScript.write("\n"+'xrdcp -s root://cmseos.fnal.gov/' + outputFolder + '/' + CMSSWRel +'.tgz  .');
outScript.write("\n"+'tar -xf '+ CMSSWRel +'.tgz' );
outScript.write("\n"+'rm '+ CMSSWRel +'.tgz' );
outScript.write("\n"+'cd ' + CMSSWRel + '/src/WWAnalysis/WWAnalysisRun2' );
outScript.write("\n"+'scramv1 b ProjectRename');
outScript.write("\n"+'eval `scram runtime -sh`');
outScript.write("\n"+command);
outScript.write("\n"+'echo "====> List output files : " ');
outScript.write("\n"+'ls WWTree*.root');
outScript.write("\n"+'xrdcp -f WWTree*.root root://cmseos.fnal.gov/' + outputFolder);
outScript.write("\n"+'rm WWTree*.root');
outScript.write("\n"+'cd ${_CONDOR_SCRATCH_DIR}');
outScript.write("\n"+'rm -rf ' + CMSSWRel);
outScript.write("\n");
outScript.close();
os.system("chmod 777 runstep2condor.sh");

outJDL = open("runstep2condor.jdl","w");
outJDL.write("Executable = runstep2condor.sh\n");
outJDL.write("Universe = vanilla\n");
#outJDL.write("Requirements =FileSystemDomain==\"fnal.gov\" && Arch==\"X86_64\"");
outJDL.write("Notification = ERROR\n");
outJDL.write("Should_Transfer_Files = YES\n");
outJDL.write("WhenToTransferOutput = ON_EXIT\n");
#outJDL.write("include : list-infiles.sh |\n");
outJDL.write("Transfer_Input_Files = "+inputlist+"\n");
outJDL.write("x509userproxy = $ENV(X509_USER_PROXY)\n");

    #MC
if( doMC ):
    for i in range(len(samples)):
        outJDL.write("Output = "+OutputLogPath+"/"+str(samples[i][1])+".stdout\n");
        outJDL.write("Error  = "+OutputLogPath+"/"+str(samples[i][1])+".stdout\n");
        outJDL.write("Arguments = -n "+str(samples[i][1])+" -o WWTree_"+str(samples[i][1])+" -w "+str(samples[i][0])+" -no "+ str(samples[i][2]) + " -noNeg " + str(samples[i][3]) + " -lumi "+str(lumi)+" --ismc 1 -trig 1 -c lpc\n");
        outJDL.write("Queue\n");
    
#data
if( doData ):
    for a in range(len(category)):
        for i in range(len(nameData[category[a]])):
            outJDL.write("Output = "+OutputLogPath+"/"+(nameData[category[a]])[i]+".stdout\n");
            outJDL.write("Error = "+OutputLogPath+"/"+(nameData[category[a]])[i]+".stdout\n");
            outJDL.write("Arguments = -n "+(nameData[category[a]])[i]+" -o WWTree_"+(nameData[category[a]])[i]+"_"+category[a]+" -w 1. -no 1. --ismc 0 -trig 1 -c lpc\n");
            outJDL.write("Queue\n");

outJDL.close();
print "===> Set Proxy Using:";
print "\tvoms-proxy-init --voms cms --valid 168:00";
print "\"condor_submit runstep2condor.jdl\" to submit";

