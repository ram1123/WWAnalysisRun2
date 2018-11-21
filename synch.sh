#!/bin/bash

v="comp"

python python/produceWWNtuples.py -i /eos/cms/store/user/jlawhorn/Synch/MC -n MC -o MC_${v} -w 1 -no 1000 --lumi 35900.0 --ismc 1 -trig 1 -c lxplus -loc 1 > MC_${v}.log
python python/produceWWNtuples.py -i /eos/cms/store/user/jlawhorn/Synch/ -n SingleEle -o SingleEle_${v} -w 1 -no 1001 --lumi 35900.0 --ismc 0 -trig 1 -c lxplus -loc 1 > SingleEle_${v}.log
python python/produceWWNtuples.py -i /eos/cms/store/user/jlawhorn/Synch/ -n SingleMu -o SingleMu_${v} -w 1 -no 1002 --lumi 35900.0 --ismc 0 -trig 1 -c lxplus -loc 1 > SingleMu_${v}.log

echo "MC---" 
diff /afs/cern.ch/work/j/jlawhorn/public/Synch/MC_init.log MC_${v}.log
echo "Electrons---" 
diff /afs/cern.ch/work/j/jlawhorn/public/Synch/SingleEle_init.log SingleEle_${v}.log
echo "Muons---"
diff /afs/cern.ch/work/j/jlawhorn/public/Synch/SingleMu_init.log SingleMu_${v}.log