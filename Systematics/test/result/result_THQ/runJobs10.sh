#!/bin/bash
export XRD_NETWORKSTACK=IPv4
export X509_USER_PROXY=/tmp/x509up_u110926
WD=$PWD
echo
echo
echo
cd /afs/cern.ch/work/p/prsaha/public/flashgg_slc7/CMSSW_9_4_9
eval $(scram runtime -sh)
cd $WD
mkdir /afs/cern.ch/work/p/prsaha/public/flashgg_slc7/CMSSW_9_4_9/src/flashgg/Systematics/test/result/result_THQ
echo "ls $X509_USER_PROXY"
ls $X509_USER_PROXY
cmsRun /afs/cern.ch/work/p/prsaha/public/flashgg_slc7/CMSSW_9_4_9/src/flashgg/Systematics/test/result/result_THQ/workspaceStd.py maxEvents=-1 campaign=RunIIFall17-3_2_0 targetLumi=41.5e+3 processType=mc processIdMap=/afs/cern.ch/work/p/prsaha/public/flashgg_slc7/CMSSW_9_4_9/src/flashgg/Systematics/test/result/result_THQ/config.json dataset=/THQ_ctcvcp_HToGG_M125_13TeV-madgraph-pythia8_TuneCP5/spigazzi-RunIIFall17-3_2_0-RunIIFall17-3_2_0-6-g0758a19d-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1-ae4aae54f34d4b504ce18c7b2946b9db/USER outputFile=/afs/cern.ch/work/p/prsaha/public/flashgg_slc7/CMSSW_9_4_9/src/flashgg/Systematics/test/result/result_THQ/output_THQ.root nJobs=38 jobId=10
retval=$?
if [[ $retval == 0 ]]; then
    errors=""
    for file in $(find -name '*.root' -or -name '*.xml'); do
        cp -pv $file /afs/cern.ch/work/p/prsaha/public/flashgg_slc7/CMSSW_9_4_9/src/flashgg/Systematics/test/result/result_THQ
        if [[ $? != 0 ]]; then
            errors="$errors $file($?)"
        fi
    done
    if [[ -n "$errors" ]]; then
       echo "Errors while staging files"
       echo "$errors"
       exit -2
    fi
fi

exit $retval

