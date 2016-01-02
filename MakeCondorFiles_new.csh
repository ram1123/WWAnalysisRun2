#!/bin/tcsh
setenv pwd $PWD
#echo ${pwd}
cat>Job_${1}.csh<<EOF
#!/bin/tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh
cd /uscms_data/d3/rasharma/WWScattering/CMSSW_7_1_5/src/WWAnalysisRun2
setenv SCRAM_ARCH slc6_amd64_gcc472
cmsenv
cd /uscms_data/d3/rasharma/WWScattering/CMSSW_7_1_5/src/WWAnalysisRun2
python python/produceWWNtuples.py -mc False -l el -n data_el_prompt_v4_25ns_runD -mass 0 -o data_el_prompt_v4_25ns_runD -trig 1
EOF

chmod 775 ${pwd}/Job_${1}.csh

cat>condor_${1}<<EOF
universe = vanilla
Executable = Job_${1}.csh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = $pwd/python/produceWWNtuples.py, $pwd/produceWWNtuples.exe, $pwd/PU.root
Output               = ${1}_\$(Cluster)_\$(Process).stdout
Error                = ${1}_\$(Cluster)_\$(Process).stderr
Log                  = ${1}_\$(Cluster)_\$(Process).log
Queue 1
EOF

condor_submit condor_${1}
condor_status rasharma
