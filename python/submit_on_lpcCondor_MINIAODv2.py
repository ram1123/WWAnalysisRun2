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
	outputFolder = "/store/user/rasharma/SecondStep/WWTree_"+datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M')+"_TEST/";
	OutputLogPath = "OutPut_Logs/Logs_" + datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M') + "_TEST";
else:
	outputFolder = "/store/user/rasharma/SecondStep/WWTree_"+datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M');
	OutputLogPath = "OutPut_Logs/Logs_" + datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M');


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
	( 6.655,	"WJetsToLNu_HT_800To1200_13TeV",	7688957,	0),
	( 6.655,	"WJetsToLNu_HT_800To1200_13TeV_ext1",	7688957,	0),
	#( 1.589,	"ZTo2LZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	99997,	0),
	#( 91.14,	"WplusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	1981107,	0),
	#( 5.633,	"WplusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	3948770,	0),
	( 0.01398,	"ZZZ_13TeV_amcatnlo_pythia8",	249232,	18020),
	#( 1.938,	"WplusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	1991348,	0),
	( 10.71,	"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_1",	23766546,	4986275),
	( 10.71,	"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_2",	23766546,	4986275),
	( 10.71,	"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_3",	23766546,	4986275),
	( 0.2529,	"TTZToLLNuNu_M-10",	7969186,	2126557),
	( 0.2529,	"TTZToLLNuNu_M-10_ext2",	7969186,	2126557),
	( 435.24,	"WJetsToLNu_HT_200To400_13TeV",	38925816,	0),
	( 435.24,	"WJetsToLNu_HT_200To400_13TeV_ext1",	38925816,	0),
	( 435.24,	"WJetsToLNu_HT_200To400_13TeV_ext2",	38925816,	0),
	( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV",	2520618,	0),
	( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV_ext1",	2520618,	0),
	( 0.2043,	"TTWJetsToLNu_ext1",	5280251,	1282079),
	( 0.2043,	"TTWJetsToLNu_ext2",	5280251,	1282079),
	( 119.9,	"QCD_HT1500to2000_13TeV",	11624720,	0),
	( 119.9,	"QCD_HT1500to2000_13TeV_ext",	11624720,	0),
	( 364.3,	"TTToSemilepton_powheg_1",	91832423,	0),
	( 364.3,	"TTToSemilepton_powheg_2",	91832423,	0),
	( 364.3,	"TTToSemilepton_powheg_3",	91832423,	0),
	( 364.3,	"TTToSemilepton_powheg_4",	91832423,	0),
	( 364.3,	"TTToSemilepton_powheg_5",	91832423,	0),
	( 364.3,	"TTToSemilepton_powheg_6",	91832423,	0),
	( 3.451,	"WplusToLNuWplusTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	199858,	0),
	#( 5.633,	"WplusTo2JWminusToLNuJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	3994663,	0),
	( 32100.0,	"QCD_HT500to700_13TeV",	53744436,	0),
	( 32100.0,	"QCD_HT500to700_13TeV_ext_1",	53744436,	0),
	( 32100.0,	"QCD_HT500to700_13TeV_ext_2",	53744436,	0),
	( 32100.0,	"QCD_HT500to700_13TeV_ext_3",	53744436,	0),
	( 27990000.0,	"QCD_HT100to200_13TeV_1",	80614044,	0),
	( 27990000.0,	"QCD_HT100to200_13TeV_2",	80614044,	0),
	( 27990000.0,	"QCD_HT100to200_13TeV_3",	80614044,	0),
	( 27990000.0,	"QCD_HT100to200_13TeV_4",	80614044,	0),
	( 27990000.0,	"QCD_HT100to200_13TeV_5",	80614044,	0),
	( 27990000.0,	"QCD_HT100to200_13TeV_6",	80614044,	0),
	( 17.94,	"WplusToLNuWminusTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	1640169,	0),
	( 3.22,	"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",	15345161,	2828391),
	( 1.895,	"WplusToLNuZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	349000,	0),
	( 5765.4,	"DYJetsToLL_M-50_amcatnlo_1",	121994032,	20127060),
	( 5765.4,	"DYJetsToLL_M-50_amcatnlo_2",	121994032,	20127060),
	( 5765.4,	"DYJetsToLL_M-50_amcatnlo_3",	121994032,	20127060),
	( 5765.4,	"DYJetsToLL_M-50_amcatnlo_4",	121994032,	20127060),
	( 5765.4,	"DYJetsToLL_M-50_amcatnlo_5",	121994032,	20127060),
	( 5765.4,	"DYJetsToLL_M-50_amcatnlo_6",	121994032,	20127060),
	( 5765.4,	"DYJetsToLL_M-50_amcatnlo_7",	121994032,	20127060),
	#( 0.03203,	"WminusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	99657,	0),
	( 11.36,	"Summer16_ST_s_channel_4f_leptonDecays",	999976,	188501),
	( 0.5686,	"WplusTo2JZTo2LJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	149363,	0),
	( 17.92,	"WplusTo2JWminusToLNuJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	1706220,	0),
	( 14.58,	"WJetsToLNu_HT_600To800_13TeV",	18578604,	0),
	( 14.58,	"WJetsToLNu_HT_600To800_13TeV_ext1",	18578604,	0),
	( 14.58,	"WJetsToLNu_HT_600To800_13TeV_ext1_1",	18578604,	0),
	( 14.58,	"WJetsToLNu_HT_600To800_13TeV_ext1_2",	18578604,	0),
	( 5.595,	"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",	26516586,	5318729),
	( 19.5741,	"ST_tW_top_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",	5372830,	0),
	( 6831.0,	"QCD_HT700to1000_13TeV",	42658677,	0),
	( 6831.0,	"QCD_HT700to1000_13TeV_ext",	42658677,	0),
	#( 3.259,	"WminusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	199525,	0),
	#( 91.07,	"WplusTo2JWminusToLNuJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	1923847,	0),
	( 1207.0,	"QCD_HT1000to1500_13TeV",	13080692,	0),
	( 1207.0,	"QCD_HT1000to1500_13TeV_ext",	13080692,	0),
	( 1627.45,	"WJetsToLNu_HT_100To200_13TeV",	79165703,	0),
	( 1627.45,	"WJetsToLNu_HT_100To200_13TeV_ext1_1",	79165703,	0),
	( 1627.45,	"WJetsToLNu_HT_100To200_13TeV_ext1_2",	79165703,	0),
	( 1627.45,	"WJetsToLNu_HT_100To200_13TeV_ext1_3",	79165703,	0),
	( 1627.45,	"WJetsToLNu_HT_100To200_13TeV_ext2_1",	79165703,	0),
	( 1627.45,	"WJetsToLNu_HT_100To200_13TeV_ext2_2",	79165703,	0),
	( 136.02,	"Summer16_ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4",	5993570,	0),
	#( 0.07584,	"WplusToLNuWplusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	99992,	0),
	( 0.5297,	"TTZToQQ",	749367,	199113),
	( 0.2223,	"WminusTo2JZTo2LJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	168679,	0),
	#( 18.25,	"WplusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	393171,	0),
	#( 0.3449,	"ZTo2LZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	49999,	0),
	( 1712000.0,	"QCD_HT200to300_13TeV",	54280031,	0),
	( 1712000.0,	"QCD_HT200to300_13TeV_ext",	54280031,	0),
	( 1712000.0,	"QCD_HT200to300_13TeV_ext_1",	54280031,	0),
	( 1712000.0,	"QCD_HT200to300_13TeV_ext_2",	54280031,	0),
	( 19.5741,	"ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",	5424956,	0),
	( 347700.0,	"QCD_HT300to500_13TeV",	26924854,	0),
	( 347700.0,	"QCD_HT300to500_13TeV_ext",	26924854,	0),
	( 0.7414,	"WminusToLNuZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	160194,	0),
	#( 0.3488,	"WminusTo2JZTo2LJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	489280,	0),
	( 59.18,	"WJetsToLNu_HT_400To600_13TeV",	7754252,	0),
	( 59.18,	"WJetsToLNu_HT_400To600_13TeV_ext1",	7754252,	0),
	#( 10.0,	"WminusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	169938,	0),
	( 10.31,	"ZZ_13TeV_pythia8",	990051,	0),
	( 0.5067,	"WminusToLNuWminusTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	190106,	0),
	#( 5.401,	"WplusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	188998,	0),
	( 1.60809,	"WJetsToLNu_HT_1200To2500_13TeV",	6708656,	0),
	( 1.60809,	"WJetsToLNu_HT_1200To2500_13TeV_ext1",	6708656,	0),
	( 49.997,	"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",	5057358,	953706),
	#( 0.575,	"WplusTo2JZTo2LJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	499432,	0),
	( 80.95,	"Summer16_ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4",	3927980,	0),
	( 0.1651,	"WWZ_13TeV_amcatnlo_pythia8",	269990,	15372),
	#( 8.793,	"WplusToLNuWplusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	198848,	0),
	( 25.24,	"QCD_HT2000toInf_13TeV",	5875869,	0),
	( 25.24,	"QCD_HT2000toInf_13TeV_ext",	5875869,	0),
	( 3.361,	"ZTo2LZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	99392,	0),
	#( 1.166,	"WminusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	981540,	0),
	#( 2.982,	"WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	198896,	0),
	( 0.4062,	"TTWJetsToQQ",	833257,	201483),
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

