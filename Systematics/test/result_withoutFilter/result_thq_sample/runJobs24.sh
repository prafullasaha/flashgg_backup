#!/bin/bash
export XRD_NETWORKSTACK=IPv4
export X509_USER_PROXY=/tmp/x509up_u110926
WD=$PWD
echo
echo
echo
cd /afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6
eval $(scram runtime -sh)
cd $WD
mkdir /afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/result_withoutFilter/result_thq_sample
echo "ls $X509_USER_PROXY"
ls $X509_USER_PROXY
cmsRun /afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/result_withoutFilter/result_thq_sample/workspaceStd.py maxEvents=-1 campaign=RunIIFall17-3_1_0 targetLumi=1e+3 processType=mc processIdMap=/afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/result_withoutFilter/result_thq_sample/config.json dataset=/THQ_ctcvcp_HToGG_M125_13TeV-madgraph-pythia8_TuneCP5/sethzenz-RunIIFall17-3_1_0-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1-042e544057588adaee2f06de8954ff07/USER outputFile=/afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/result_withoutFilter/result_thq_sample/output_THQ.root nJobs=38 jobId=24
retval=$?
if [[ $retval == 0 ]]; then
    errors=""
    for file in $(find -name '*.root' -or -name '*.xml'); do
        cp -pv $file /afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/result_withoutFilter/result_thq_sample
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

