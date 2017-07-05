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
    parser.add_argument ('-i', '--inputFolder' , default = '/store/user/arapyan/Run2/WWJJToLNuQQ_LL_13TeV-madgraph-pythia8/Samples/170503_175516/' , help='input folder with the reduced trees')
    #parser.add_argument ('-i', '--inputFolder' , default = '/store/user/arapyan/Run2/WWJJToLNuQQ_LT_13TeV-madgraph-pythia8/Samples/170503_175158/' , help='input folder with the reduced trees')
    #parser.add_argument ('-i', '--inputFolder' , default = '/store/cmst3/group/monojet/production/' , help='input folder with the reduced trees')
    parser.add_argument ('-o', '--output' , default = 'OutPutRootFile', help='output file')
    parser.add_argument ('-v', '--vbfsel' , default = '2', help='1 = select highest pt jet pair, 2 = select highest mjj, 3 = select highest DEta_jj VBF Jets')
    parser.add_argument ('-mc', '--ismc' , default = '0', help='is MC or not')
    parser.add_argument ('-l', '--lepton' , default = 'el', help='lepton category (mu or el)')
    parser.add_argument ('-t', '--tree' , default = 'Events', help='name of the input tree')
    parser.add_argument ('-n', '--name' , default = '0000' , help='input file')
    parser.add_argument ('-w', '--xsecWeight' , default = '0.0002739' , help='xsec (pb)')
    parser.add_argument ('-no', '--numberOfEntries' , default = '28687' , help='number of initial entries of the dataset')
    parser.add_argument ('-mass', '--mass' , default = '1000' , help='mass of input signal')
    parser.add_argument ('-trig', '--applyTrigger' , default = '0' , help='apply trigger or not')
    parser.add_argument ('-json', '--json', default = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt', help="json file to apply")
    parser.add_argument ('-loc', '--isLocal' , default = '0', help='run in local or not')
    parser.add_argument ('-exe', '--exe' , default = './produceWWNtuples.exe', help='location of the executable')
    args = parser.parse_args ()

    command = args.exe+' '+args.inputFolder+' '+args.output+' '+args.ismc+' '+args.lepton+' '+args.tree+' '+args.name+' '+args.xsecWeight+' '+args.numberOfEntries+' '+args.mass+' '+args.applyTrigger+' '+args.json+' '+args.isLocal+' '+args.vbfsel
    print command
    os.system(command)
