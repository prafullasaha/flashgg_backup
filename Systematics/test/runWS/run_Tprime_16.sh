fggRunJobs.py --load $CMSSW_BASE/src/flashgg/Systematics/test/samples/Legacy16/sample16_Tprime.json -d sig_16 --stage-to /eos/user/p/prsaha/Tprime_analysis/ws/output_Tprime_16 -x cmsRun ../workspaceStd.py -n 50 -q testmatch --no-copy-proxy useAAA=1 dumpTrees=false dumpWorkspace=true
#copyInputMicroAOD=1
