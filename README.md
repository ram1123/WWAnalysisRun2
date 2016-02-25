# WWAnalysisRun2

The package contains a code to produce ntuple for WW semileptonic final state.
It takes in input ntuples produced from miniAOD with the TreeMaker at this link: https://github.com/lbrianza/RA2_2014


Instructions:

git clone https://github.com/lbrianza/WWAnalysisRun2;

Few points to note:
	1. The main code is in directory bin named produceWWNtuples.cpp
	2. This takes a input text file that should be in pwd.
		1. This input text file has list of input root file (with path) for a particular sample.
		2. Presently you can find all the input text file in directory InputRootFiles

cd WWAnalysisRun2/;

make;

python python/produceWWNtuples.py -n ReducedSelection_RSGraviton4000.root -o RSGraviton4000.root -mc True
