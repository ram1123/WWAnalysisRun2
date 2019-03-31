#
#	Before running this script need to modify two things:
#	
#	1. Copy file runstep2condor.jdl to runstep2condor.jdl.bk
#	2. Name of files in Argument line
#	3. Run this script. But the ResubTemplate1 part will be at the end of file. So, cut and paste it to top of file.
#

ResubTemplate1='''
Executable = runstep2condor.sh
Universe = vanilla
Notification = ERROR
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = runstep2condor.sh, python/produceWWNtuples.py
x509userproxy = $ENV(X509_USER_PROXY)
'''
ResubTemplate2='''
Output = OutPut_Logs/Logs_2017_12_19_12h35/%s.stdout
Error  = OutPut_Logs/Logs_2017_12_19_12h35/%s.stdout
Arguments = -n %s -o WWTree_%s -w %s -no %s -noNeg %s -lumi 35900.0 --ismc 1 -trig 1 -c lpc
Queue
'''

Argument = [
"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_8",
"DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_3",
"DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_5",
"DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1",
"DYToLL_2J_13TeV-amcatnloFXFX-pythia8_1",
"WplusTo2JWminusToLNuJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",
"WplusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",
"WplusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",
"WplusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",
"WplusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8"
]


print ResubTemplate1,

import os
for args in Argument:
	os.system('grep '+args+' runstep2condor.jdl.bk')
	os.system('echo Queue')
