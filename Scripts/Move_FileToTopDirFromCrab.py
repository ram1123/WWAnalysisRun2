#!/usr/bin/python

import sys
import os

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 2:
        print "please enter path of input directory as one of argumetns of command\n"
	print "python Move_FileToTopDirFromCrab.py  <Path_of_Top_Crab_Directory>\n"
        exit(0)


source = sys.argv[1]


#os.system("eos root://cmseos.fnal.gov ls "+sys.argv[1])

search1='log'
search2='failed'
Arrayfilepath = []
for root, dirs, filenames in os.walk(source):   
        for f in filenames:
                filepath = root + os.sep + f
                if filepath.find(search2) == -1:        # Don't select file if it contains word failed in path
                        if filepath.find(search1) == -1:        # Don't select file if it contains wored log in path
                                if filepath.endswith(".root"):  # select only root files
                                        Arrayfilepath.append(filepath)
#print len(Arrayfilepath)                                       
#print Arrayfilepath[0]

for i in range(0,len(Arrayfilepath)):
        print "====> Copying file:",i," => ",Arrayfilepath[i]
        #print("xrdcp -f "+ Arrayfilepath[i]+" "+ source)
        os.system("xrdcp -f "+ Arrayfilepath[i]+" "+ source)
