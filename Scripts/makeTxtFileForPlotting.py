#!/usr/bin/python
import os
import re
import yaml
import timeit
import uproot 
from colorpicker import * 
from hadd_files_from_array import hadd_root_files


import argparse

PARSER = argparse.ArgumentParser(description='User inputs')

PARSER.add_argument('-d', '--hadd',
		    action='store_false', 
		    default=False,
		    help='Tell if you need to perform only hadd on the similar root files or not (default: False, it means it will perform hadd as well as it will generate the text file.)'
		   ) # ifhaddOnly
PARSER.add_argument('-p', '--path',
		    action='store', 
		    default='/store/user/rasharma/SecondStep/Run_2017/Frameworkupdate/WWTree_2019_05_01_14h50/', 
		    required=True,
		    help='Path of store are where it will find the input root file (default = /store/user/rasharma/SecondStep/Run_2017/Frameworkupdate/WWTree_2019_05_01_14h50/).'
		   ) # StoreArea
PARSER.add_argument('-s', '--save',
		    action='store', 
		    default='/uscms_data/d3/rasharma/aQGC_analysis/PlottingMacros/CMSSW_9_0_1/src/PlottingCodes2017/ControlPlots/SampleFiles/',
		    help='Add path of plotting code direcotory. After making the text file it will copy the same to this directory (default = /uscms_data/d3/rasharma/aQGC_analysis/PlottingMacros/CMSSW_9_0_1/src/PlottingCodes2017/ControlPlots/SampleFiles/).'
		   )  # PlottingDirectoryPath
PARSER.add_argument('-a', '--arrange',
		    action='store', 
		    default=['# name', 'data', 'Data', 'WV_EWK', 'aQGC', 'Diboson', 'VV', 'W\+jets', 'Z\+jets', 'top', 'QCD'], 
		    help="Arrange the each line of the generated text file in the given order (default = ['# name', 'data', 'Data', 'WV_EWK', 'aQGC', 'Diboson', 'VV', 'W\+jets', 'Z\+jets', 'top', 'QCD'])"
		   ) # searchString
PARSER.add_argument('-o', '--outdir',
		    action='store', 
		    default='HaddedFiles',
		    help='Name of output directory where it will place the hadd-ed root files (default = HaddedFiles)'
		   ) # StoreAreaHadd = StoreArea+'HaddedFiles_Test/'
PARSER.add_argument('-e', '--eosstring',
		    action='store', 
		    default='/eos/uscms',
		    help='Add the store area initials using which one can access it locally (default = /eos/uscms).'
		   ) #source = "/eos/uscms"+StoreArea
ARGS = PARSER.parse_args()

#print ARGS

def check_directory(path, directory_name):
    """Check if output directory exists or not. 
       - If the output directory does not exists, create it.
       - If it exists then remove the output directory, then
         create the new one.
    """
    print "path = ", path
    print "directory name = ", directory_name
    if os.path.isdir(ARGS.eosstring+path+os.sep+directory_name):
        print(CRED+"Directory "+path+os.sep+directory_name+' found. Delete it...'+CEND)
        os.system('eos root://cmseos.fnal.gov/ rm -r '+path+os.sep+directory_name)
        print "\n# list all directories: Just to check if directory HaddedFiles directory does not exists"
        os.system('ls -ltrh '+ARGS.eosstring+path+' | grep ^d')
        print "\n\nCreate the directory: ", path+os.sep+directory_name
        os.system('eos root://cmseos.fnal.gov/ mkdir '+path+os.sep+directory_name)
    else:
        print "Create the directory"
        os.system('eos root://cmseos.fnal.gov/ mkdir '+path+os.sep+directory_name)
    print "\n", "#"*51, "\n\n"

def getfilename_array(path, fileextension):
    """Get an array containing all root files
    from the "path"
    """
    arrayfile_path = []
    count = 0
    print "\nList all root file in ", ARGS.eosstring+path
    for root, dirs, filenames in os.walk(ARGS.eosstring+path):
        for filepath in sorted(filenames):
            if filepath.endswith(fileextension):
                arrayfile_path.append(filepath)
                count += 1
                print count, filepath
    return sorted(arrayfile_path)

def loadyaml_file(filename):
    """Load the yaml file info as a dict.
    """
    infile = open(filename,"r")
    ymload = yaml.load(infile)
    infile.close()
    return ymload

def getfiles_attributes(yamlfilecontent, filename_withpath_array, path):
    """Check sample information in Yaml file and 
    search corresponding sample name in the path
    then count number of total events present in the 
    path and the number of negative events
    """
    #TODO: Reverse the looping condition and if any root file
    #	    present at path is not found in the Yaml file then
    #	    report the error
    root_file_lists_array = []
    event_count_list = []
    count_negative_events_list = []
    for sample in sorted(yamlfilecontent):
        temp = []
        temp.append(sample+".root")
        count_events = 0
        count_negative_events = 0
        for i, files in enumerate(filename_withpath_array):
            if files.find(sample) != -1:
                print CBOLD+"====================> "+sample+CEND
                print CGREEN, "File-Name\t\tnEvents\tnNegative-Events", CEND
                #print i, files
                if os.path.isfile(ARGS.eosstring+path+os.sep+files):
                    print CRED, i, files, CEND
		    #if files.find('Single') != -1:
                    #FIXME: change the reading method to follow eos guideline
                    if (uproot.open(ARGS.eosstring+path+os.sep+files).keys()) == []:
                        print (CGREEN+"\nskip file: "+files+"\n"+CEND)
                    else:
                        #FIXME: change the reading method from eos area
                        otree = uproot.open(ARGS.eosstring+path+os.sep+files)["otree"]
                        inputarrays = otree.arrays(["nEvents", "nNegEvents"])
                        if len(inputarrays["nEvents"]):
                            #TODO: Here we can improve speed for data, by skipping earlier if there is data. As we don't count number of events.
                            if files.find('Single') != -1:
                                count_events = 1
                                count_negative_events = 0
                            else:
                                if ARGS.hadd != 1:
                                    count_events += inputarrays["nEvents"][0]
                                    count_negative_events += inputarrays["nNegEvents"][0] 
                                else:
                                    print "skip..."
                            temp.append(files)
                            print files, count_events, count_negative_events
                else:
                    print(CRED+"\nFile Not Found: "+files+"\n"+CEND)
        if len(temp)>1:
            root_file_lists_array.append(temp)
            event_count_list.append(count_events)
            count_negative_events_list.append(count_negative_events)
    return root_file_lists_array, event_count_list, count_negative_events_list

def haddfiles_listoflists(input_dir, output_dir, listoflists_rootfiles):
    """Perform the hadd operation for the list of input files.
    """
    for root_file_array in listoflists_rootfiles:
        temp = "hadd -f "
        for count_file in range(0, len(root_file_array)):
            if count_file == 0:
                temp += output_dir+os.sep+root_file_array[count_file]+' '
            else:
                temp += input_dir+os.sep+root_file_array[count_file]+" "
        print temp
        os.system(temp)

def create_plotting_textfile(filepath, root_file_lists_array, ymload,
                             event_count_list, count_negative_events_list,
			     outputfile):
    """Generate a text file have summarized info of file name, its
    number of events, number of negative events, cross-section.
    This file we need as an input to the plotting macro
    """
    print "=============  MAKE SUMMARY    ================"  
    print(root_file_lists_array)
    
    print outputfile
    write_textfile = open(outputfile, "w")
    write_textfile.write("# name           file_location  xspb/lumipb  otherscale nMCevents       nMCNegEvents    colorcode       stackit\n")
    print ("# name           file_location  xspb/lumipb  otherscale nMCevents       nMCNegEvents    colorcode       stackit\n")
    for i, files in enumerate(root_file_lists_array):
        print i, files[0]
        for sample in ymload:
            if files[0].find(sample) != -1:
                print i, files[0]
                print ymload[sample]["name"], "\t", filepath+os.sep+files[0], "\t", ymload[sample]["CrossSection"], "\t1\t", int(event_count_list[i]), "\t", int(count_negative_events_list[i]), "\t", ymload[sample]["ColorCode"], "\t", ymload[sample]["StackIt"]
                write_textfile.write(ymload[sample]["name"]+"\t"+filepath+os.sep+files[0]+"\t"+str(ymload[sample]["CrossSection"])+"\t1\t"+str(event_count_list[i])+"\t"+str(count_negative_events_list[i])+"\t"+str(ymload[sample]["ColorCode"])+"\t"+str(ymload[sample]["StackIt"])+"\n")
    write_textfile.close()

def sorted_textfile(outputfile):
    """Sort the text file as required by the plotting macro
    """
    sorted_outputfile = open("Sorted_"+outputfile,'w')
    
    for strings in ARGS.arrange:
        infile = open(outputfile, "r")
        for lines in infile:
            if re.match(str(strings), lines): 
                sorted_outputfile.write(lines)
    sorted_outputfile.close()
  
################################################

def main():
    """Main program for generating the summarized text file
    for the plotting macro.
    """
    print CBOLD+CGREEN+"Step: 1: Check if the input and output directory exists"+CEND+CEND
    check_directory(ARGS.path, ARGS.outdir)

    print CBOLD+CGREEN+"Step: 2: Make list of all root files in eos."+CEND+CEND
    filename_withpath_array = getfilename_array(ARGS.path, ".root")

    print CBOLD+CGREEN+"Step: 3: Read Yaml file."+CEND+CEND
    yamlfilecontent = loadyaml_file("DataMCInfo.yml")

    print CBOLD+CGREEN+"\n\nStep: 4: Read each root file and saves its number of events and number of negative events in each file."+CEND+CEND
    print CBOLD+CGREEN+"Step: 5: Club the root files based on the keys given in the yaml file"+CEND+CEND
    root_file_lists_array, count_events, count_negative_events = getfiles_attributes(yamlfilecontent, filename_withpath_array, ARGS.path)
    print root_file_lists_array, count_events, count_negative_events

    if ARGS.hadd:
        print CBOLD+CGREEN+"\n\nStep: 6: Do the hadd for each clubbed list."+CEND+CEND
        #haddfiles_listoflists(ARGS.eosstring+ARGS.path, ARGS.eosstring+ARGS.path+os.sep+ARGS.outdir, root_file_lists_array)
	hadd_root_files(root_file_lists_array, ARGS.eosstring+ARGS.path, ARGS.eosstring+ARGS.path+os.sep+ARGS.outdir)
    else:
        print CBOLD+CRED+"\n\nDon't perform the hadd step."+CEND+CEND

    print CBOLD+CGREEN+"\n\nStep: 7: Create the text file that contains all info for the plotting."+CEND+CEND
    outputfile = "DibosonBoostedElMuSamples13TeV_"+(ARGS.eosstring+ARGS.path+os.sep+ARGS.outdir).replace("/","_")+".txt"
    create_plotting_textfile(ARGS.path+os.sep+ARGS.outdir, root_file_lists_array, yamlfilecontent, count_events, count_negative_events, outputfile)

    print CBOLD+CGREEN+"\n\nStep: 8: Sort the output text file."+CEND+CEND
    sorted_textfile(outputfile)

    # copy the created text file to the respective directory
    print CBOLD+CGREEN+"\n\nStep: 9: Copy the plotting text file to output directory: ", ARGS.path+os.sep+ARGS.outdir+CEND+CEND
    os.system('xrdcp -f '+outputfile+'  root://cmseos.fnal.gov/'+ ARGS.path+os.sep+ARGS.outdir)
    os.system('xrdcp -f '+"Sorted_"+outputfile+'  root://cmseos.fnal.gov/'+ ARGS.path+os.sep+ARGS.outdir)


if __name__ == "__main__":
    START_CLOCK = timeit.default_timer() # start timer
    main()
    STOP_CLOCK = timeit.default_timer() # stop timer
    # Print total time to run in minutes
    print 'Run Time: ', round((STOP_CLOCK - START_CLOCK)/60.0, 2), 'min'  

