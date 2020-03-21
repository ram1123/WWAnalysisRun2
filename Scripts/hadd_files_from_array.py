import os
import commands
import argparse
import timeit

def hadd_files_from_dir(input_path, output_path, search_string):
    print input_path
    print output_path
    print search_string
    print "#"*51
    rootFileList = []
    rootFileList.append("output-temp.root")
    for root, dirs, filenames in os.walk("/eos/uscms"+input_path):
	for f in filenames:
	    if search_string in f: 
	        #print f
		rootFileList.append(f)
    print "size of list = ",len(rootFileList)
    hadd_root_files(rootFileList, input_path, output_path)

def hadd_root_files(input_root_files, input_path, output_path):
    print "Number of root files = ",len(input_root_files)
    simple_hadd(input_root_files, input_path, output_path)
    """
    if len(input_root_files) < ARGS.nRootfiles:
        simple_hadd(input_root_files, input_path, output_path)
    else:
        #hadd_many_files(input_root_files, input_path, output_path)
        simple_hadd(input_root_files, input_path, output_path)
    """

def simple_hadd(input_root_files, input_path, output_path):
    #print "Input files is < "+ARGS.nRootfiles+", so, performing simple hadd"
    print "="*51
    print "\n\n"
    print "="*51
    print "input_root_files[0] = \n",input_root_files[0]
    print "="*51
    print "input_path = ",input_path
    print "="*51
    print "output_path = ",output_path
    print "="*51
    #command = "hadd -f root://cmsxrootd.fnal.gov/"+output_path.replace("/eos/uscms","")+os.sep+input_root_files[0]
    command = "hadd -f "+output_path+os.sep+input_root_files[0]
    for root_file in range(1,len(input_root_files)):
        command += " "+input_path+os.sep+input_root_files[root_file]
        #command += " "+"/eos/uscms/"+input_path+os.sep+input_root_files[root_file]
    print command
    os.system(command)
    print "done."

def hadd_many_files(input_root_files, input_path, output_path):
    hadd_command = ""
    tempCount = 0
    os.system("rm -r /uscmst1b_scratch/lpc1/3DayLifetime/"+commands.getstatusoutput('echo ${USER}')[1]+"/tmp")
    os.system("mkdir -p /uscmst1b_scratch/lpc1/3DayLifetime/"+commands.getstatusoutput('echo ${USER}')[1]+"/tmp/haddfiles")
    temp_directory = "/uscmst1b_scratch/lpc1/3DayLifetime/"+commands.getstatusoutput('echo ${USER}')[1]+"/tmp/haddfiles"
    print "Temporary directory for placing root files: ",temp_directory
    os.system("mkdir -p "+temp_directory)
    for count_files, root_file in enumerate(input_root_files):
        if count_files == 0:
	    pass
	else:
	    if count_files%ARGS.nRootfiles == 1:
	        print "\n\n","#"*51
		if tempCount != 0:
		    print "Hadd command:\n\t",hadd_command
		    os.system(hadd_command)
		else:
		    print "Don't print; Reading first file\n"
		hadd_command = "hadd -f "+temp_directory+os.sep+"temp_"+str(tempCount)+"_temp.root "+"/eos/uscms"+input_path+os.sep+root_file
		tempCount += 1
	    else:
	        hadd_command = hadd_command + " " + "/eos/uscms"+input_path+os.sep+root_file
    print "\n\n","#"*51
    print "Hadd command:\n\t",hadd_command
    os.system(hadd_command)
    
    hadd_command=""
    for root, dirs, filenames in os.walk(temp_directory):
        rootFileList = []
	for f in filenames:
	    if f.endswith("_temp.root"): rootFileList.append(f)
	    tempCount=0
	    print rootFileList
	    hadd_command = "hadd -f "+"/eos/uscms"+output_path+os.sep+input_root_files[0]
            for count,rootFile in enumerate(rootFileList):
                print "\n\n","#"*51
                hadd_command = hadd_command + " " + temp_directory+os.sep+rootFile
        del rootFileList[:]
    print hadd_command
    os.system(hadd_command)
    print "done..."

def main():
    """hadd files
    """
    #print "input array: ", ARGS.input_array
    print "\n\npath: ",ARGS.path
    print "\n\noutput path: ",ARGS.output_path
    #hadd_root_files(ARGS.input_array, ARGS.path, ARGS.output_path)
    hadd_files_from_dir(ARGS.path, ARGS.output_path, ARGS.string_to_search)

if __name__ == "__main__":
    """Hadd files
    """
    PARENT_PARSER = argparse.ArgumentParser(add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    
    PARENT_PARSER.add_argument('-p', '--path',
    		    action='store', 
    		    required=True,
    		    help='Path of store are where it will find the input root file'
    		   ) # StoreArea
    PARENT_PARSER.add_argument('-o', '--output_path',
    		    action='store', 
    		    required=True,
    		    help='Path of store are where it will find the input root file'
    		   ) # StoreAre
    PARENT_PARSER.add_argument('-n', '--nRootfiles',
	  action='store',
	  default=108,
	  help='number of root files to hadd at a time'
	  )
    
    PARSER = argparse.ArgumentParser(description='User inputs')
    
    SUBPARSER = PARSER.add_subparsers(help="Available subcommands and their descriptions.")
    
    PARSE_HADD_FILES = SUBPARSER.add_parser("hadd_files", help='', parents=[PARENT_PARSER])
    
    PARSE_HADD_DIR = SUBPARSER.add_parser("hadd_dir", help='', parents=[PARENT_PARSER], formatter_class=argparse.RawTextHelpFormatter)
    
    PARSE_HADD_FILES.add_argument('-i', '--input_array',
    		    action='store', 
    		    required=True,
    		    nargs='+',
    		    #type=str,
    		    #dest='list',
    		    default=[],
    		    help="Input array should be of type:['out.root', 'in_1.root', 'in_2.root', 'in_3.root']"
    		   ) # ifhaddOnly
    PARSE_HADD_DIR.add_argument('-s', '--string_to_search',
                        action='store',
    		    required=True,
    		    default="test",
    		    help="string to search"
    		   )
    
    ARGS = PARSER.parse_args()

    print "#"*51
    START_CLOCK = timeit.default_timer() # start timer
    main()
    STOP_CLOCK = timeit.default_timer() # stop timer
    # Print total time to run in minutes
    print 'Run Time: ', round((STOP_CLOCK - START_CLOCK)/60.0, 2), 'min'  

