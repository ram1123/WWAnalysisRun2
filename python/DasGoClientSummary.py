import os
from subprocess import Popen, PIPE

datasets = [	"/*EWK_LO_SM*/*TrancheIV*/MINIAODSIM",
		"/*QCD_LO_SM*/*TrancheIV*/MINIAODSIM",
		"/*EWK_LO_aQGC*/*TrancheIV*/MINIAODSIM",
		"/DYToLL_0J_13TeV-amcatnloFXFX-pythia8*/*TrancheIV*/MINIAODSIM",
		"/DYToLL_1J_13TeV-amcatnloFXFX-pythia8*/*TrancheIV*/MINIAODSIM",
		"/DYToLL_2J_13TeV-amcatnloFXFX-pythia8*/*TrancheIV*/MINIAODSIM",
		"/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*/*TrancheIV*/MINIAODSIM",
		"/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*/*TrancheIV*/MINIAODSIM",
		"/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*/*TrancheIV*/MINIAODSIM",
		"/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*/*TrancheIV*/MINIAODSIM",
		#***"/DYJetsToLL_M-50_T*amcatnlo*/*/MINIAODSIM",
		"/DYJetsToLL_M-50_T*amcatnlo*pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X*mcRun2*TrancheIV*/MINIAODSIM",
		"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8*/*TrancheIV*/MINIAODSIM",
		"/WJetsToLNu_HT_100To200*13TeV*/*TrancheIV*/MINIAODSIM",
		"/WJetsToLNu_HT_200To400*13TeV*/*TrancheIV*/MINIAODSIM",
		"/WJetsToLNu_HT_400To600*13TeV*/*TrancheIV*/MINIAODSIM",
		"/WJetsToLNu_HT_600To800*13TeV*/*TrancheIV*/MINIAODSIM",
		"/WJetsToLNu_HT_800To1200*13TeV*/*TrancheIV*/MINIAODSIM",
		"/WJetsToLNu_HT_1200To2500*13TeV*/*TrancheIV*/MINIAODSIM",
		"/WJetsToLNu_HT_2500ToInf*13TeV*/*TrancheIV*/MINIAODSIM",
		"/WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8*/*TrancheIV*/MINIAODSIM",
		#**"/ZZ_13TeV_pythia8*/*/MINIAODSIM",
		"/ZZZ_TuneCUETP8M1_13TeV*pythia8*/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
		"/ZZ_TuneCUETP8M1_13TeV*pythia8*/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
		"/ZZTo2L2Q_13TeV_powheg_pythia8*/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
		"/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8*/*TrancheIV*/MINIAODSIM",
		"/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8*/*TrancheIV*/MINIAODSIM",
		"/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8*/*TrancheIV*/MINIAODSIM",
		"/WWZ*13TeV*amcatnlo*pythia8*/*TrancheIV*/MINIAODSIM",
		"/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
		"/TTZToLLNuNu_M-10*/*TrancheIV*/MINIAODSIM",
		"/TTZToQQ*/*TrancheIV*/MINIAODSIM",
		"/TTZToLLNuNu_M-10*/*TrancheIV*/MINIAODSIM",
		"/TTWJetsToQQ*/*TrancheIV*/MINIAODSIM",
		"/TTWJetsToLNu*/*TrancheIV*/MINIAODSIM",
		"/ST_s_channel*4f*leptonDecays*/RunIISummer16MiniAODv2*TrancheIV*/MINIAODSIM",
		"/ST_tW_top_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1*/*TrancheIV*/MINIAODSIM",
		"/ST_t-channel*inclusiveDecays*TuneCUETP8M2T4*/RunIISummer16MiniAODv2*TrancheIV*/MINIAODSIM",
		"/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1*/*TrancheIV*/MINIAODSIM",
		"/QCD_HT*TuneCUETP8M1_13TeV*/RunIISummer16MiniAODv2*TrancheIV*/MINIAODSIM"
		#**"/QCD_HT100to200_13TeV*/*/MINIAODSIM",
		#**"/QCD_HT200to300_13TeV*/*/MINIAODSIM",
		#**"/QCD_HT300to500_13TeV*/*/MINIAODSIM",
		#**"/QCD_HT500to700_13TeV*/*/MINIAODSIM",
		#**"/QCD_HT700to1000_13TeV*/*/MINIAODSIM",
		#**"/QCD_HT1000to1500_13TeV*/*/MINIAODSIM",
		#**"/QCD_HT1500to2000_13TeV*/*/MINIAODSIM",
		#**"/QCD_HT2000toInf_13TeV*/*/MINIAODSIM"
		#"/Single*/*3Feb2017*/MINIAOD"
	] 

count1=1
for i,datas in enumerate(datasets):
	#command1 = 'dasgoclient --query="dataset='+datas+'" --limit=0'
	#print command1
	#print ""
	pipe = Popen('dasgoclient --query="dataset='+datas+'" --limit=0', shell=True, stdout=PIPE)
	#pipe = Popen(command1, shell=True, stdout=PIPE)

	count2=1
	for line in pipe.stdout:
		print "Job Nubmer = ",count1,":",count2
    		print "Dataset: ",line.strip(),'\n'
    		command = 'dasgoclient --query="summary dataset='+line.strip()+'"'
    		os.system(command)
    		print "\n-----------\n"
		count2+=1
	count1+=1
