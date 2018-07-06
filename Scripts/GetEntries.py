import os
import ROOT as ROOT

ROOT.gROOT.SetBatch(True)

source = '/uscms_data/d3/rasharma/aQGC_analysis/SecondStep/CMSSW_8_0_26_patch1/src/WWAnalysis/WWAnalysisRun2/Output_4Sep/'

Arrayfilepath = []
for root, dirs, filenames in os.walk(source):	
	for f in filenames:
		#filepath = root + os.sep + f
		filepath =  f
		if filepath.endswith(".root"):
			Arrayfilepath.append(filepath)

for i in range(0,len(Arrayfilepath)):
	#file = ROOT.TFile(Arrayfilepath[i])
	file = ROOT.TFile(source+Arrayfilepath[i])
	#print Arrayfilepath[i]
	if not file:
		print '\n==>Failed to open %s' % Arrayfilepath[i]
		print '\n'
		exit(0)
	tree = file.Get('otree')
	h1 = ROOT.TH1F("h1","", 999999999, 0 , 999999999)
	h2 = ROOT.TH1F("h2","", 999999999, 0 , 999999999)
	tree.Draw("nNegEvents>>h1","","",10)
	tree.Draw("nEvents>>h2","","",10)
	#print Arrayfilepath[i],"\t",tree.GetEntries()
	print Arrayfilepath[i],"\t",int(h2.GetMean()),"\t",int(h1.GetMean())
	h1.Delete()
	h2.Delete()
