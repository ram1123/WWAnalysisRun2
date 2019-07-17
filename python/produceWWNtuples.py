#!/usr/bin/python
import sys
import os
import commands
from commands import getstatusoutput
import datetime
import argparse
import string
if __name__ == '__main__':
    parser = argparse.ArgumentParser (description = 'produce ntuples with WW semileptonic final state')
    parser.add_argument ('-i', '--inputFolder' , default = '/store/user/arapyan/Run2/' , help='input folder with the reduced trees')
    parser.add_argument ('-o', '--output' , default = 'OutPutRootFile', help='output file')
    parser.add_argument ('-v', '--vbfsel' , default = '2', help='1 = select highest pt jet pair, 2 = select highest mjj, 3 = select highest DEta_jj VBF Jets')
    parser.add_argument ('-mc', '--ismc' , default = '0', help='is MC or not')
    parser.add_argument ('-c', '--cluster' , default = 'lxplus', help='cluster can be lxplus or lpc')
    parser.add_argument ('-t', '--tree' , default = 'Events', help='name of the input tree')
    parser.add_argument ('-n', '--name' , default = 'WplusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8' , help='input file')
    parser.add_argument ('-f', '--infile' , default = 'ttHTobb_M125_TuneCP5_13TeV_powheg_pythia8_0.txt' , help='input text file having path of root files to run over')
    parser.add_argument ('-w', '--xsecWeight' , default = '0.0002739' , help='xsec (pb)')
    parser.add_argument ('-no', '--numberOfEntries' , default = '28687' , help='number of initial entries of the dataset')
    parser.add_argument ('-noNeg', '--numberOfNegEntries' , default = '0' , help='number of initial entries of the negative events in dataset')
    parser.add_argument ('-lumi', '--lumi' , default = '35900' , help='Luminosity in pb-1')
    parser.add_argument ('-trig', '--applyTrigger' , default = '0' , help='apply trigger or not')
    parser.add_argument ('-json', '--json', default = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt', help="json file to apply")
    parser.add_argument ('-loc', '--isLocal' , default = '0', help='run in local or not')
    parser.add_argument ('-exe', '--exe' , default = 'produceWWNtuples', help='location of the executable')
    args = parser.parse_args ()

    command = args.exe+' '+args.inputFolder+' '+args.output+' '+args.ismc+' '+args.cluster+' '+args.tree+' '+args.infile+' '+args.xsecWeight+' '+args.numberOfEntries+' '+args.numberOfNegEntries+' '+args.applyTrigger+' '+args.json+' '+args.isLocal+' '+args.vbfsel
    print "#"*51,"\n#\n#\tCommand to run\n#\n# \t",
    print command,
    print "\n#\n","#"*51
    os.system(command)
