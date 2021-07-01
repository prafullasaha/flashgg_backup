fggRunJobs.py --load $CMSSW_BASE/src/flashgg/Systematics/test/samples/Legacy17/sample_data17.json -d data_17 --stage-to /eos/user/p/prsaha/Tprime_analysis/ws/output_data_17 -x cmsRun ../workspaceStd.py -n 50 -q testmatch --no-copy-proxy dumpTrees=false dumpWorkspace=true
#copyInputMicroAOD=1
