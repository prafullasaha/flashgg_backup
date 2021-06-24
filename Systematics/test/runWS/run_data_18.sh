fggRunJobs.py --load $CMSSW_BASE/src/flashgg/Systematics/test/samples/Legacy18/sample_data18.json -d data_18 --stage-to /eos/user/p/prsaha/Tprime_analysis/ws/output_data_18 -x cmsRun ../workspaceStd.py -n 50 -q testmatch --no-copy-proxy dumpTrees=false dumpWorkspace=true
#copyInputMicroAOD=1
