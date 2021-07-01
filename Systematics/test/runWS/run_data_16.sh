fggRunJobs.py --load $CMSSW_BASE/src/flashgg/Systematics/test/samples/Legacy16/sample_data16.json -d data_16 --stage-to /eos/user/p/prsaha/Tprime_analysis/ws/output_data_16 -x cmsRun ../workspaceStd.py -n 50 -q testmatch --no-copy-proxy dumpTrees=false dumpWorkspace=true
#copyInputMicroAOD=1
