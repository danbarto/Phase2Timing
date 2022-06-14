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
git remote add upstream git@github.com:danbarto/Phase2Timing.git
```

Compile
```
scram build clean
scram b -j 8
```

