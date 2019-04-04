#!/bin/bash

#source setup.sh
#kinit
for x in `grep -H "Fatal in <mithep::Selector::AbortAnalysis>: Analysis aborted!\|Error in <TXNetFile::CreateXClient>: open attempt failed on root\|segmentation violation\|Job Killed.\|Unable to get quota space - quota not defined or exhausted\|bus error\|Exited with exit code 255.\|Exited\|is not defined in current scope" out.* | sed "s@:@ @g" | awk '{print $1}' | uniq |  head -500`;  do 
#for x in `grep -H -L "Info in <mithep::Analysis::Run>: Processing complete" out.* | sed "s@:@ @g" | awk '{print $1}' | uniq |  head -500`;  do 
    line=`grep 'Job <' $x | awk  '{print $2 }' | sed 's/<//g' | sed 's/>//g'`
    #key=`echo $line | awk '{print $3}'`
    #newkey=`ls -ltr /tmp/ | grep krb | grep \`id -u\` | sed "s@_@ @g" | tail -1 | awk '{print $11}'`
    #line=`echo $line | sed "s@$key@$newkey@g"`
    bsub -o out.%J -q 2nd $line
    rm $x
    echo $line
    #folder=`echo $x | sed "s@\/@ @g" | awk '{print $1 }'`
    #echo $folder
    #rm -r $folder
done

