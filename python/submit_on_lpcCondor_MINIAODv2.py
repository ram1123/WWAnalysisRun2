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

inputFolder = "/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/";
TestRun = 0

doMC = True;
doData = True;
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


print "Name of output dir: ",outputFolder
# create a directory on eos
os.system('xrdfs root://cmseos.fnal.gov/ mkdir ' + outputFolder)
# create directory in pwd for log files
os.system('mkdir -p ' + OutputLogPath)

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
os.system('xrdcp -f ' + CMSSWRel+".tgz" + ' root://cmseos.fnal.gov//store/user/rasharma/' + CMSSWRel+".tgz")
#with open("ThingsUpdated.txt","a") as myfile:
#	CmdOutput = subprocess.check_output('git log --pretty=format:"%h - %an, %cd : %s" -2')
#	myfile.write("=========================================\n\n")
#	myfile.write(CmdOutput)
#myfile.close()
os.system('xrdcp -f ThingsUpdated.txt root://cmseos.fnal.gov/' + outputFolder)
os.system('cp ThingsUpdated.txt ' + OutputLogPath)

samples = [
    ##( 0.9114,	"WplusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ##( 0.9107,	"WplusTo2JWminusToLNuJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ##( 0.0879,	"WplusToLNuWplusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ##( 0.0326,	"WminusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ##( 0.1825,	"WplusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
    ##( 0.0540,	"WplusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
    ##( 0.1000,	"WminusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
    ##( 0.0298,	"WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
    ##( 0.0159,	"ZTo2LZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
    ##( 5.546,	"WplusTo2JWminusToLNuJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ##( 5.558,	"WplusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ##( 0.086,	"WplusToLNuWplusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ##( 0.038,	"WminusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ##( 2.159,	"WplusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
    ##( 0.640,	"WplusTo2JZTo2LJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
    ##( 1.302,	"WminusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
    ##( 0.387,	"WminusTo2JZTo2LJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
    ##( 0.376,	"ZTo2LZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
    ( 17.94,	"WplusToLNuWminusTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ( 17.92,	"WplusTo2JWminusToLNuJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ( 3.451,	"WplusToLNuWplusTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ( 0.5067,	"WminusToLNuWminusTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ( 1.895,	"WplusToLNuZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ( 0.5686,	"WplusTo2JZTo2LJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
    ( 0.7414,	"WminusToLNuZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ( 0.2223,	"WminusTo2JZTo2LJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ( 3.361,	"ZTo2LZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
    ( 542.00,	"TTToSemilepton_powheg_1",	0),
    ( 542.00,	"TTToSemilepton_powheg_2",	0),
    ( 542.00,	"TTToSemilepton_powheg_3",	0),
    ( 542.00,	"TTToSemilepton_powheg_4",	0),
    ( 542.00,	"TTToSemilepton_powheg_5",	0),
    ( 542.00,	"TTToSemilepton_powheg_6",	0),
    #( ,	"TTTo2L2Nu_powheg",	0),
    #( ,	"TTTo2L2Nu_powheg_pythia8",	0),
    ( 0.2043,	"TTWJetsToLNu_ext1",	0),
    ( 0.2043,	"TTWJetsToLNu_ext2",	0),
    ( 0.4062,	"TTWJetsToQQ",	0),
    ( 0.2529,	"TTZToLLNuNu_M-10",	0),
    ( 0.2529,	"TTZToLLNuNu_M-10_ext2",	0),
    ( 0.5297,	"TTZToQQ",	0),
    #( ,	"ST_tW_antitop_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4",	0),
    #( ,	"ST_tW_top_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4",	0),
    ( 19.47,	"ST_tW_top_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",	0),
    ( 19.47,	"ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",	0),
    ( 11.36,	"Summer16_ST_s_channel_4f_leptonDecays",	0),
    ( 80.95,	"Summer16_ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4",	0),
    ( 136.02,	"Summer16_ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4",	0),
    #( ,	"WW_13TeV_pythia8",	0),
    ( 10.31,	"ZZ_13TeV_pythia8",	0),
    ( 0.01398,	"ZZZ_13TeV_amcatnlo_pythia8",	0),
    ( 0.1651,	"WWZ_13TeV_amcatnlo_pythia8",	0),
    ( 49.997,	"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",	0),
    ( 3.22,	"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",	0),
    ( 5.595,	"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",	0),
    ( 10.71,	"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_1",	0),
    ( 10.71,	"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_2",	0),
    ( 10.71,	"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_3",	0),
    ( 49.997,	"WWToLNuQQ_13TeV_powheg",			0),
    ( 49.997,	"WWToLNuQQ_13TeV_powheg_ext",			0),
    #( ,	"WWW_4F_13TeV_amcatnlo_pythia8",	0),
    #( ,	"ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8",	0),
    #( ,	"ZZto4L_13TeV_amcatnloFXFX_pythia8",	0),
    #( ,	"ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8",	0),
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo_1",	0),
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo_2",	0),
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo_3",	0),
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo_4",	0),
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo_5",	0),
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo_6",	0),
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo_7",	0),
    #( 1.68,	"DYJetsToLL_M-50_HT-600to800",	0),
    ( 1378.3414,	"WJetsToLNu_HT_70To100_13TeV",	0),
    ( 1506.4,		"WJetsToLNu_HT_100To200_13TeV", 		0),
    ( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext1_1",		0),
    ( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext1_2",		0),
    ( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext1_3",		0),
    ( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext2_1",		0),
    ( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext2_2",		0),
    ( 435.237,		"WJetsToLNu_HT_200To400_13TeV",			0),
    ( 435.237,		"WJetsToLNu_HT_200To400_13TeV_ext1",		0),
    ( 435.237,		"WJetsToLNu_HT_200To400_13TeV_ext2",		0),
    ( 59.1811,		"WJetsToLNu_HT_400To600_13TeV",			0),
    ( 59.1811,		"WJetsToLNu_HT_400To600_13TeV_ext1",		0),
    ( 14.5805,		"WJetsToLNu_HT_600To800_13TeV",			0),
    ( 14.5805,		"WJetsToLNu_HT_600To800_13TeV_ext1_1",		0),
    ( 14.5805,		"WJetsToLNu_HT_600To800_13TeV_ext1_2",		0),
    ( 6.65621,		"WJetsToLNu_HT_800To1200_13TeV",		0),
    ( 6.65621,		"WJetsToLNu_HT_800To1200_13TeV_ext1",		0),
    ( 1.60809,		"WJetsToLNu_HT_1200To2500_13TeV",		0),
    ( 1.60809,		"WJetsToLNu_HT_1200To2500_13TeV_ext1",		0),
    ( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV",		0),
    ( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV_ext1",		0),
    ( 256380272,	"QCD_HT50to100_13TeV",	0),
    ( 27990000,		"QCD_HT100to200_13TeV_1",		0),
    ( 27990000,		"QCD_HT100to200_13TeV_2",		0),
    ( 27990000,		"QCD_HT100to200_13TeV_3",		0),
    ( 27990000,		"QCD_HT100to200_13TeV_4",		0),
    ( 27990000,		"QCD_HT100to200_13TeV_5",		0),
    ( 27990000,		"QCD_HT100to200_13TeV_6",		0),
    ( 1712000,		"QCD_HT200to300_13TeV",		0),
    ( 1712000,		"QCD_HT200to300_13TeV_ext_1",	0),
    ( 1712000,		"QCD_HT200to300_13TeV_ext_2",	0),
    ( 347700,		"QCD_HT300to500_13TeV",		0),
    ( 347700,		"QCD_HT300to500_13TeV_ext",	0),
    ( 32100,		"QCD_HT500to700_13TeV",		0),
    ( 32100,		"QCD_HT500to700_13TeV_ext_1",	0),
    ( 32100,		"QCD_HT500to700_13TeV_ext_2",	0),
    ( 32100,		"QCD_HT500to700_13TeV_ext_3",	0),
    ( 6831,		"QCD_HT700to1000_13TeV",	0),
    ( 6831,		"QCD_HT700to1000_13TeV_ext",	0),
    ( 1207,		"QCD_HT1000to1500_13TeV",	0),
    ( 1207,		"QCD_HT1000to1500_13TeV_ext",	0),
    ( 119.9,		"QCD_HT1500to2000_13TeV",	0),
    ( 119.9,		"QCD_HT1500to2000_13TeV_ext",	0),
    ( 25.24,		"QCD_HT2000toInf_13TeV",	0),
    ( 25.24,		"QCD_HT2000toInf_13TeV_ext",	0),
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
outScript.write("\n"+'xrdcp -s root://cmseos.fnal.gov//store/user/rasharma/' + CMSSWRel +'.tgz  .');
outScript.write("\n"+'tar -xf '+ CMSSWRel +'.tgz' );
outScript.write("\n"+'rm '+ CMSSWRel +'.tgz' );
outScript.write("\n"+'cd ' + CMSSWRel + '/src/WWAnalysis/WWAnalysisRun2' );
outScript.write("\n"+'scramv1 b ProjectRename');
outScript.write("\n"+'eval `scram runtime -sh`');
outScript.write("\n"+command);
outScript.write("\n"+'echo "====> List output files : " ');
outScript.write("\n"+'ls WWTree*.root');
outScript.write("\n"+'xrdcp WWTree*.root root://cmseos.fnal.gov/' + outputFolder);
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
        outJDL.write("Arguments = -n "+str(samples[i][1])+" -o WWTree_"+str(samples[i][1])+" -w "+str(samples[i][0])+" -lumi "+str(lumi)+" --ismc 1 -trig 1 -c lpc\n");
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

