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
mkdir /afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/result_cutflow/with_pho_pt
echo "ls $X509_USER_PROXY"
ls $X509_USER_PROXY
cmsRun /afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/result_cutflow/with_pho_pt/workspaceStd.py maxEvents=-1 campaign=RunIIFall17-3_1_0 targetLumi=1e+3 processType=mc processIdMap=/afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/result_cutflow/with_pho_pt/config.json dataset=/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/sethzenz-RunIIFall17-3_1_0-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1-6a727f257a984f9b93ec84d5e8acd977/USER outputFile=/afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/result_cutflow/with_pho_pt/output_TTHJet.root nJobs=10 jobId=6
retval=$?
if [[ $retval == 0 ]]; then
    errors=""
    for file in $(find -name '*.root' -or -name '*.xml'); do
        cp -pv $file /afs/cern.ch/work/p/prsaha/public/mytest1/CMSSW_9_4_6/src/flashgg/Systematics/test/result_cutflow/with_pho_pt
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

