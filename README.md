# Phase2Timing
Repo for testing phase 2 timing methods

See https://twiki.cern.ch/twiki/bin/view/CMS/FLT for details

## Setup

Setup CMSSW:
```
cmsrel CMSSW_11_3_1_patch1 
cd CMSSW_11_3_1_patch1/src/
cmsenv
```

Fork this repo and then
```
git clone git@github.com:<YOUR_USER>/Phase2Timing.git
cd Phase2Timing
git remote add upstream git@github.com:danbarto/Phase2Timing.git
```

Compile
```
scram build clean
scram b -j 8
```

## Run locally
Don't forget about `cmsenv`!

Inside the `Phase2Timing` directory do
```
cmsRun Phase2TimingAnalyzer/python/ConfFile_local_file_cfg.py
```

Input signal files are on ceph, e.g.:

``` shell
/ceph/cms//store/user/mcitron/ProjectMetis/HTo2LongLivedTo4e_MH-125_MFF-50_CTau-1000mm_privateMC_11X_RECOMINI_v1_generationForPhase2HS_noPU_CEPH_vector/
/ceph/cms//store/user/mcitron/ProjectMetis/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-1000mm_privateMC_11X_RECOMINI_v1_generationForPhase2_noPU_CEPH/
```


