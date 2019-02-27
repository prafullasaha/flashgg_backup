#!/usr/bin/env sh

cd ${_CONDOR_SCRATCH_DIR}
echo ${_CONDOR_SCRATCH_DIR}
source /cvmfs/cms.cern.ch/cmsset_default.sh
#echo "xrdcp root://eoscms.cern.ch//afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/CMSSW_9_4_6.tar ."
cmsrel CMSSW_9_4_6
cd CMSSW_9_4_6/src
scp /afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/CMSSW_9_4_6_src.tar.gz .
echo "Copied"
echo "tar -xzf CMSSW_9_4_6_src.tar.gz"
tar -xzf CMSSW_9_4_6_src.tar.gz
echo "file extracted"
echo $ls -lhtr
#scp /afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/external .
#source /cvmfs/cms.cern.ch/cmsset_default.sh

echo "inside src directory"
eval `scram runtime -sh`
scram build -j 8
echo "Compilation complete"
source flashgg/afterbuild_9_4_X.sh
cd flashgg/Systematics/test/
echo "Inside test directory"
#sleep 5
#eval `scramv1 runtime -sh`

echo "fggRunJobs.py --load sample_thq.json -D -P -n 100 -d root://eoscms.cern.ch//afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/result_fortalk/result_thq -x cmsRun workspaceStd.py maxEvents=-1 "

fggRunJobs.py --load sample_thq.json -D -P -n 100 -d /afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/result_fortalk/result_thq -x cmsRun workspaceStd.py maxEvents=-1

